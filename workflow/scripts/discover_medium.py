"""Discover and replay one MAG's minimum-total-import growth medium."""

import json
import math
import sys
import traceback
from pathlib import Path

sys.path.insert(0, str(snakemake.scriptdir))

from common.media import cobra_nutrient_mappings, write_medium_csv
from cplex_setup import ensure_cplex_ready


def _log(message):
    text = str(message)
    print(text, flush=True)
    log_path = Path(snakemake.log[0])
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a", encoding="utf-8") as handle:
        handle.write(text + "\n")


def _set_threads(model, threads):
    problem = getattr(model.solver, "problem", None)
    parameters = getattr(problem, "parameters", None)
    if parameters is not None and hasattr(parameters, "threads"):
        parameters.threads.set(int(threads))


def main():
    model_path = str(snakemake.input.model)
    output = str(snakemake.output.medium)
    summary_out = str(snakemake.output.summary)
    min_growth = float(snakemake.params.min_growth)
    open_exchanges = bool(snakemake.params.open_exchanges)
    exchange_bound = float(snakemake.params.exchange_bound)
    solver = str(snakemake.params.solver)
    mock = bool(snakemake.params.mock)
    mag = str(snakemake.wildcards.mag)

    with Path(snakemake.input.filtered_models).open(encoding="utf-8") as handle:
        retained_mags = {line.strip() for line in handle if line.strip()}
    if mag not in retained_mags:
        raise ValueError(
            f"MAG {mag!r} is not present in the post-MEMOTE retained cohort "
            f"declared by {snakemake.input.filtered_models}."
        )

    if mock:
        medium = {"mock_nutrient_e": 1.0}
        write_medium_csv(output, medium)
        Path(summary_out).parent.mkdir(parents=True, exist_ok=True)
        Path(summary_out).write_text(
            json.dumps(
                {
                    "synthetic": True,
                    "mag": mag,
                    "source": "discover",
                    "min_growth": min_growth,
                    "open_exchanges": open_exchanges,
                    "exchange_bound": exchange_bound,
                    "verified_growth": min_growth,
                    "nutrient_count": 1,
                },
                indent=2,
                sort_keys=True,
            )
            + "\n",
            encoding="utf-8",
        )
        _log("Synthetic mock medium written; no optimization was performed.")
        return

    if solver != "cplex":
        raise ValueError(
            "Native simulation-medium discovery is supported with solver=cplex "
            "in the pinned MICOM/COBRApy environment."
        )
    ensure_cplex_ready(
        solver,
        snakemake.config.get("cplex", {}).get("studio_dir", ""),
        auto_upgrade=bool(snakemake.config.get("cplex", {}).get("auto_upgrade", True)),
        log_file=str(snakemake.log[0]),
    )

    from cobra.io import read_sbml_model
    from cobra.medium import minimal_medium

    model = read_sbml_model(model_path)
    model.solver = solver
    _set_threads(model, snakemake.threads)
    _log(
        f"Discovering minimum-total-import medium for MAG={mag} "
        f"min_growth={min_growth} open_exchanges={open_exchanges} "
        f"exchange_bound={exchange_bound}."
    )
    recipe = minimal_medium(
        model,
        min_objective_value=min_growth,
        exports=False,
        minimize_components=False,
        open_exchanges=exchange_bound if open_exchanges else False,
    )
    if recipe is None or recipe.empty:
        raise RuntimeError(
            f"No minimum-import medium supports growth {min_growth} for MAG "
            f"{mag!r}; check its objective, exchange reactions, and configured "
            "discovery opening settings."
        )

    all_exchange_mappings = cobra_nutrient_mappings(
        model,
        model_path,
        [
            metabolite.id
            for reaction in model.exchanges
            for metabolite in reaction.metabolites
        ],
    )
    reaction_to_mapping = {}
    for mapping in all_exchange_mappings.values():
        if mapping.cobra_exchange_id in reaction_to_mapping:
            raise ValueError(
                f"Exchange {mapping.cobra_exchange_id!r} in {model_path} maps to "
                "multiple nutrient identities."
            )
        reaction_to_mapping[mapping.cobra_exchange_id] = mapping

    medium = {}
    replay = {}
    for exchange_id, uptake in recipe.items():
        uptake = float(uptake)
        if not math.isfinite(uptake) or uptake <= 0:
            raise RuntimeError(
                f"Discovery returned invalid uptake {uptake!r} for exchange "
                f"{exchange_id!r} in MAG {mag!r}."
            )
        if exchange_id not in reaction_to_mapping:
            raise ValueError(
                f"Discovery exchange {exchange_id!r} in MAG {mag!r} cannot be "
                "mapped to one actual extracellular nutrient identity."
            )
        mapping = reaction_to_mapping[exchange_id]
        if mapping.nutrient_id in medium:
            raise ValueError(
                f"MAG {mag!r} has ambiguous duplicate import exchanges for "
                f"nutrient {mapping.nutrient_id!r}."
            )
        medium[mapping.nutrient_id] = uptake
        replay[exchange_id] = uptake

    if open_exchanges:
        # minimal_medium applies this opening only inside its temporary model
        # context. Recreate the same candidate-exchange assumptions before the
        # exact recipe is applied so replay does not fall back to narrower
        # source secretion bounds. The medium setter then closes every
        # nonselected uptake while retaining the configured secretion side.
        for reaction in model.exchanges:
            reaction.bounds = (-exchange_bound, exchange_bound)
    model.medium = replay
    solution = model.optimize()
    tolerance = float(model.solver.configuration.tolerances.feasibility)
    growth = float(solution.objective_value)
    if solution.status != "optimal" or not math.isfinite(growth):
        raise RuntimeError(
            f"Replaying the discovered medium for MAG {mag!r} did not produce "
            f"a finite optimal solution (status={solution.status!r}, "
            f"growth={growth!r})."
        )
    if growth + tolerance < min_growth:
        raise RuntimeError(
            f"Replaying the discovered medium for MAG {mag!r} achieved growth "
            f"{growth}, below requested {min_growth} with solver tolerance "
            f"{tolerance}."
        )

    write_medium_csv(output, medium)
    Path(summary_out).parent.mkdir(parents=True, exist_ok=True)
    Path(summary_out).write_text(
        json.dumps(
            {
                "synthetic": False,
                "mag": mag,
                "source": "discover",
                "solver": solver,
                "threads": int(snakemake.threads),
                "min_growth": min_growth,
                "open_exchanges": open_exchanges,
                "exchange_bound": exchange_bound,
                "verification_status": str(solution.status),
                "verified_growth": growth,
                "solver_tolerance": tolerance,
                "nutrient_count": len(medium),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    _log(f"Published {len(medium)} nutrient(s); replay growth={growth:.12g}.")


try:
    main()
except Exception:
    _log("Medium discovery failed:\n" + traceback.format_exc())
    raise
