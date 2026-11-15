"""Run native MICOM FBA and pFBA on the one selected simulation medium."""

import json
import math
import sys
import traceback
from pathlib import Path

import pandas as pd

SCRIPT_DIR = (
    Path(snakemake.scriptdir)
    if "snakemake" in globals()
    else Path(__file__).resolve().parent
)
sys.path.insert(0, str(SCRIPT_DIR))

from common.abundance import parse_and_normalize_abundances
from common.media import (
    backend_ids,
    model_ids_and_paths,
    read_mapping_tsv,
    read_medium_csv,
)
from cplex_setup import ensure_cplex_ready


def _log(message):
    text = str(message)
    print(text, flush=True)
    log_path = Path(snakemake.log[0])
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a", encoding="utf-8") as handle:
        handle.write(text + "\n")


def _one_path(value):
    if isinstance(value, (str, Path)):
        return str(value)
    values = list(value)
    return str(values[0]) if values else ""


def _write_tsv(path, rows, columns, *, na_rep=""):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows, columns=columns).to_csv(
        path, sep="\t", index=False, na_rep=na_rep
    )


def _abundances(model_ids):
    abundance_path = _one_path(snakemake.input.get("abundance", ""))
    if abundance_path:
        normalized = parse_and_normalize_abundances(abundance_path, model_ids)
        if normalized.extra_ids:
            _log(
                "Ignoring abundance IDs without Memote-passing models: "
                f"{list(normalized.extra_ids)}"
            )
        values = list(normalized.values)
        table_format = normalized.table_format
    else:
        values = [1.0 / len(model_ids)] * len(model_ids)
        table_format = "equal"
        _log("No MICOM abundance TSV provided; assigning equal abundances.")
    zero_ids = [model_id for model_id, value in zip(model_ids, values) if value == 0.0]
    positive_ids = [
        model_id for model_id, value in zip(model_ids, values) if value > 0.0
    ]
    if zero_ids:
        _log(
            "Retaining zero-abundance MAGs without pseudocounts: "
            f"{zero_ids}; positive-abundance MAGs={positive_ids}."
        )
    return values, table_format, zero_ids


def _validate_status(path, tag):
    status = json.loads(Path(path).read_text(encoding="utf-8"))
    if status.get("tag") != tag:
        raise ValueError(
            f"Selected-medium status tag {status.get('tag')!r} does not match "
            f"requested MICOM tag {tag!r}."
        )
    return status


def _set_threads(community, threads):
    problem = getattr(community.solver, "problem", None)
    parameters = getattr(problem, "parameters", None)
    if parameters is not None and hasattr(parameters, "threads"):
        parameters.threads.set(int(threads))


def _micom_medium(community, mapping_rows, medium, raw_to_backend, abundance_by_id):
    """Resolve canonical nutrients through actual member and medium objects."""
    pool_by_nutrient = {}
    zero_unlinked = []
    for row in mapping_rows:
        raw_id = row["mag"]
        if raw_id not in raw_to_backend:
            raise ValueError(
                f"Medium mapping references MAG {raw_id!r} outside the retained cohort."
            )
        nutrient = row["nutrient_id"]
        if nutrient not in medium:
            raise ValueError(
                f"Medium mapping references unselected nutrient {nutrient!r}."
            )
        backend = raw_to_backend[raw_id]
        reactions = [
            reaction
            for reaction in community.reactions
            if getattr(reaction, "community_id", None) == backend
            and getattr(reaction, "global_id", None) == row["cobra_exchange_id"]
        ]
        if len(reactions) != 1:
            raise ValueError(
                f"MICOM member {raw_id!r} exposes {len(reactions)} actual "
                f"exchange objects for {row['cobra_exchange_id']!r}; expected one."
            )
        reaction = reactions[0]
        member_metabolites = [
            metabolite
            for metabolite in reaction.metabolites
            if getattr(metabolite, "community_id", None) == backend
            and getattr(metabolite, "global_id", None) == nutrient
        ]
        if len(member_metabolites) != 1:
            raise ValueError(
                f"MICOM member exchange {reaction.id!r} does not map uniquely to "
                f"actual nutrient object {nutrient!r}."
            )
        coefficient = float(reaction.metabolites[member_metabolites[0]])
        expected = float(row["exchange_coefficient"])
        if not math.isclose(coefficient, expected, rel_tol=0.0, abs_tol=1e-12):
            raise ValueError(
                f"MICOM changed source exchange scaling for {raw_id!r}/"
                f"{nutrient!r}: source={expected}, member={coefficient}."
            )

        medium_metabolites = [
            metabolite
            for metabolite, value in reaction.metabolites.items()
            if getattr(metabolite, "community_id", None) == "medium"
            and abs(float(value)) > 1e-15
        ]
        if not medium_metabolites and abundance_by_id[raw_id] == 0.0:
            zero_unlinked.append((raw_id, nutrient))
            continue
        if len(medium_metabolites) != 1:
            raise ValueError(
                f"MICOM member exchange {reaction.id!r} maps to "
                f"{len(medium_metabolites)} actual medium objects; expected one."
            )
        medium_metabolite = medium_metabolites[0]
        outer = [
            exchange
            for exchange in community.exchanges
            if getattr(exchange, "community_id", None) == "medium"
            and medium_metabolite in exchange.metabolites
        ]
        if len(outer) != 1:
            raise ValueError(
                f"MICOM nutrient {nutrient!r} maps to {len(outer)} actual outer "
                "exchange objects; expected one."
            )
        previous = pool_by_nutrient.setdefault(nutrient, outer[0])
        if previous is not outer[0]:
            raise ValueError(
                f"MICOM nutrient {nutrient!r} collides across multiple medium pools."
            )

    reverse = {}
    for nutrient, exchange in pool_by_nutrient.items():
        previous = reverse.setdefault(exchange.id, nutrient)
        if previous != nutrient:
            raise ValueError(
                f"MICOM outer exchange {exchange.id!r} collides between nutrients "
                f"{previous!r} and {nutrient!r}."
            )
    if not pool_by_nutrient:
        raise ValueError(
            "No selected nutrient maps through an actual MICOM member/medium "
            "exchange path."
        )
    if zero_unlinked:
        _log(
            "Zero-abundance member exchanges have no medium coefficient and were "
            f"not used to infer pools: {zero_unlinked}."
        )
    selected = {
        exchange.id: medium[nutrient] for nutrient, exchange in pool_by_nutrient.items()
    }
    community.medium = selected
    active = community.medium
    unexpected = sorted(set(active) - set(selected))
    missing = sorted(set(selected) - set(active))
    if unexpected or missing:
        raise RuntimeError(
            "MICOM did not apply the selected medium exactly "
            f"(unexpected active imports={unexpected}, missing imports={missing})."
        )
    return {nutrient: exchange.id for nutrient, exchange in pool_by_nutrient.items()}


def _validate_solution_metrics(stage, status, growth, objective):
    status = str(status)
    growth = float(growth)
    objective = float(objective)
    if status != "optimal" or not all(map(math.isfinite, (growth, objective))):
        raise RuntimeError(
            f"MICOM {stage} failed to return a finite optimal solution "
            f"(status={status!r}, growth={growth!r}, "
            f"objective={objective!r})."
        )
    return growth, objective


def _solution_rows(stage, solution, backend_to_raw):
    growth, objective = _validate_solution_metrics(
        stage,
        solution.status,
        solution.growth_rate,
        solution.objective_value,
    )
    expected = set(backend_to_raw)
    available = set(solution.members.index) & expected
    if available != expected:
        raise RuntimeError(
            f"MICOM {stage} member results do not contain the complete cohort "
            f"(missing={sorted(expected - available)})."
        )

    community_row = {
        "stage": stage,
        "growth": growth,
        "objective_value": objective,
        "stage_status": "complete",
        "synthetic": False,
    }
    member_rows = []
    for backend, raw_id in backend_to_raw.items():
        row = solution.members.loc[backend]
        abundance = float(row["abundance"])
        member_growth = float(row["growth_rate"])
        if not math.isfinite(abundance) or abundance < 0:
            raise RuntimeError(
                f"MICOM {stage} returned invalid abundance {abundance!r} for "
                f"MAG {raw_id!r}."
            )
        if not math.isfinite(member_growth):
            raise RuntimeError(
                f"MICOM {stage} returned invalid growth {member_growth!r} for "
                f"MAG {raw_id!r}."
            )
        member_rows.append(
            {
                "stage": stage,
                "mag": raw_id,
                "abundance": abundance,
                "growth": member_growth,
                "synthetic": False,
            }
        )

    flux_rows = []
    for compartment, row in solution.fluxes.iterrows():
        output_compartment = backend_to_raw.get(compartment, str(compartment))
        for reaction_id, value in row.items():
            if pd.isna(value):
                continue
            flux = float(value)
            if not math.isfinite(flux):
                raise RuntimeError(
                    f"MICOM {stage} returned non-finite flux for "
                    f"{compartment!r}/{reaction_id!r}."
                )
            flux_rows.append(
                {
                    "stage": stage,
                    "compartment": output_compartment,
                    "reaction_id": reaction_id,
                    "flux": flux,
                    "synthetic": False,
                }
            )
    if not flux_rows:
        raise RuntimeError(f"MICOM {stage} returned no available flux values.")
    return community_row, member_rows, flux_rows


def _solver_diagnostic(community, solution=None):
    if solution is None:
        return {"native_status": None, "native_message": None}
    status = getattr(solution, "status", None)
    message = None
    problem = getattr(getattr(community, "solver", None), "problem", None)
    native_solution = getattr(problem, "solution", None)
    get_status_string = getattr(native_solution, "get_status_string", None)
    if callable(get_status_string):
        try:
            message = str(get_status_string())
        except Exception:
            message = None
    if message == "Unknown status value":
        message = None
    return {
        "native_status": None if status is None else str(status),
        "native_message": message,
    }


def _complete_stage_status(diagnostic):
    return {
        "result_status": "complete",
        **diagnostic,
        "exception_type": None,
        "exception_message": None,
    }


def _failed_stage_status(diagnostic, error):
    return {
        "result_status": "unavailable",
        **diagnostic,
        "exception_type": type(error).__name__,
        "exception_message": str(error),
    }


def _pfba_failure_warning(stage_status):
    return (
        "WARNING: PARTIAL MICOM results: validated FBA outputs are retained, but "
        "supplementary pFBA is unavailable "
        f"(native_status={stage_status['native_status']!r}, "
        f"native_message={stage_status['native_message']!r}, "
        f"exception={stage_status['exception_type']}: "
        f"{stage_status['exception_message']}). Failed pFBA numerical values are "
        "not exported."
    )


def _run_micom_analyses(community, backend_to_raw, log=_log):
    fba_solution = community.optimize(fluxes=True, pfba=False, raise_error=True)
    fba_row, fba_members, fba_fluxes = _solution_rows(
        "fba", fba_solution, backend_to_raw
    )
    stage_status = {
        "fba": _complete_stage_status(_solver_diagnostic(community, fba_solution))
    }
    community_rows = [fba_row]
    member_rows = list(fba_members)
    flux_rows = list(fba_fluxes)

    pfba_solution = None
    try:
        pfba_solution = community.optimize(fluxes=True, pfba=True, raise_error=True)
        pfba_row, pfba_members, pfba_fluxes = _solution_rows(
            "pfba", pfba_solution, backend_to_raw
        )
    except Exception as error:
        diagnostic = _solver_diagnostic(community, pfba_solution)
        stage_status["pfba"] = _failed_stage_status(diagnostic, error)
        community_rows.append(
            {
                "stage": "pfba",
                "growth": None,
                "objective_value": None,
                "stage_status": "unavailable",
                "synthetic": False,
            }
        )
        log(_pfba_failure_warning(stage_status["pfba"]))
        analysis_status = "partial"
    else:
        stage_status["pfba"] = _complete_stage_status(
            _solver_diagnostic(community, pfba_solution)
        )
        community_rows.append(pfba_row)
        member_rows.extend(pfba_members)
        flux_rows.extend(pfba_fluxes)
        analysis_status = "complete"
    return (
        community_rows,
        member_rows,
        flux_rows,
        analysis_status,
        stage_status,
    )


def _publish_results(output, community_rows, member_rows, flux_rows, status_payload):
    _write_tsv(
        output.community_growth,
        community_rows,
        ["stage", "growth", "objective_value", "stage_status", "synthetic"],
        na_rep="NA",
    )
    _write_tsv(
        output.member_growth,
        member_rows,
        ["stage", "mag", "abundance", "growth", "synthetic"],
    )
    _write_tsv(
        output.fluxes,
        flux_rows,
        ["stage", "compartment", "reaction_id", "flux", "synthetic"],
    )
    Path(output.status).parent.mkdir(parents=True, exist_ok=True)
    Path(output.status).write_text(
        json.dumps(status_payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def main():
    tag = str(snakemake.params.tag)
    mock = bool(snakemake.params.mock)
    solver = str(snakemake.params.solver)
    model_ids, model_paths = model_ids_and_paths(
        snakemake.input.models, snakemake.input.filtered_models
    )
    values, abundance_format, zero_ids = _abundances(model_ids)
    medium = read_medium_csv(snakemake.input.medium)
    medium_status = _validate_status(snakemake.input.medium_status, tag)
    raw_to_backend = backend_ids(model_ids)
    backend_to_raw = {backend: raw for raw, backend in raw_to_backend.items()}

    if mock:
        community_rows = []
        member_rows = []
        flux_rows = []
        for stage in ("fba", "pfba"):
            community_rows.append(
                {
                    "stage": stage,
                    "growth": 0.0,
                    "objective_value": 0.0,
                    "stage_status": "complete",
                    "synthetic": True,
                }
            )
            member_rows.extend(
                {
                    "stage": stage,
                    "mag": model_id,
                    "abundance": abundance,
                    "growth": 0.0,
                    "synthetic": True,
                }
                for model_id, abundance in zip(model_ids, values)
            )
            flux_rows.extend(
                {
                    "stage": stage,
                    "compartment": "medium",
                    "reaction_id": f"synthetic_import_{nutrient}",
                    "flux": -uptake,
                    "synthetic": True,
                }
                for nutrient, uptake in medium.items()
            )
        pool_map = {nutrient: f"synthetic_import_{nutrient}" for nutrient in medium}
        analysis_status = "complete"
        stage_status = {
            stage: {
                "result_status": "complete",
                "native_status": None,
                "native_message": None,
                "exception_type": None,
                "exception_message": None,
            }
            for stage in ("fba", "pfba")
        }
    else:
        if solver != "cplex":
            raise ValueError(
                "Native MICOM simulation is supported with solver=cplex in the "
                "pinned MICOM/COBRApy environment."
            )
        ensure_cplex_ready(
            solver,
            snakemake.config.get("cplex", {}).get("studio_dir", ""),
            auto_upgrade=bool(
                snakemake.config.get("cplex", {}).get("auto_upgrade", True)
            ),
            log_file=str(snakemake.log[0]),
        )
        from micom import Community

        taxonomy = pd.DataFrame(
            {
                "id": [raw_to_backend[model_id] for model_id in model_ids],
                "file": model_paths,
                "abundance": values,
            }
        )
        community = Community(
            taxonomy,
            solver=solver,
            progress=False,
            rel_threshold=-1.0,
            max_exchange=1000,
        )
        if set(community.taxa) != set(backend_to_raw):
            raise RuntimeError(
                "MICOM did not retain the complete model cohort; no models may be "
                "silently dropped."
            )
        _set_threads(community, snakemake.threads)
        pool_map = _micom_medium(
            community,
            read_mapping_tsv(snakemake.input.medium_mapping),
            medium,
            raw_to_backend,
            dict(zip(model_ids, values)),
        )
        _log(
            f"Running MICOM FBA and pFBA for {len(model_ids)} model(s) on "
            f"selected medium {tag!r}."
        )
        (
            community_rows,
            member_rows,
            flux_rows,
            analysis_status,
            stage_status,
        ) = _run_micom_analyses(community, backend_to_raw)

    status_payload = {
        "synthetic": mock,
        "tool": "micom",
        "medium_tag": tag,
        "medium_source": medium_status.get("source"),
        "solver": solver,
        "threads": int(snakemake.threads),
        "analysis_status": analysis_status,
        "stages": ["fba", "pfba"],
        "stage_status": stage_status,
        "model_ids": model_ids,
        "abundance_format": abundance_format,
        "zero_abundance_model_ids": zero_ids,
        "medium_exchange_by_nutrient": pool_map,
    }
    _publish_results(
        snakemake.output,
        community_rows,
        member_rows,
        flux_rows,
        status_payload,
    )


if "snakemake" in globals():
    try:
        main()
    except Exception:
        _log("MICOM failed:\n" + traceback.format_exc())
        raise
