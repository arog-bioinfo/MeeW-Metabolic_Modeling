"""Run native ReFramed SteadyCom and SteadyComVA on the selected medium."""

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


def _write_tsv(path, rows, columns, *, na_rep=""):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows, columns=columns).to_csv(
        path, sep="\t", index=False, na_rep=na_rep
    )


def _validate_status(path, tag):
    status = json.loads(Path(path).read_text(encoding="utf-8"))
    if status.get("tag") != tag:
        raise ValueError(
            f"Selected-medium status tag {status.get('tag')!r} does not match "
            f"requested SteadyCom tag {tag!r}."
        )
    return status


def _medium_constraints(community, mapping_rows, medium, raw_to_backend):
    merged = community.merged_model
    exchange_ids = set(merged.get_exchange_reactions())
    environmental_by_nutrient = {}
    raw_to_nutrient = {}
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
        raw_species = row["raw_species_id"]
        previous = raw_to_nutrient.setdefault(raw_species, nutrient)
        if previous != nutrient:
            raise ValueError(
                f"Raw species {raw_species!r} collides between nutrients "
                f"{previous!r} and {nutrient!r}."
            )
        organism = community.organisms[raw_to_backend[raw_id]]
        if raw_species not in organism.metabolites:
            raise ValueError(
                f"ReFramed model for MAG {raw_id!r} does not contain mapped raw "
                f"species {raw_species!r}."
            )
        if raw_species not in organism.get_external_metabolites():
            raise ValueError(
                f"Mapped raw species {raw_species!r} for MAG {raw_id!r} is not "
                "extracellular in the loaded ReFramed model."
            )
        if raw_species not in merged.metabolites:
            raise ValueError(
                f"Merged ReFramed community lacks actual environmental species "
                f"{raw_species!r}."
            )
        candidates = [
            reaction_id
            for reaction_id in merged.get_metabolite_reactions(raw_species)
            if reaction_id in exchange_ids
        ]
        if len(candidates) != 1:
            raise ValueError(
                f"Merged ReFramed species {raw_species!r} maps to "
                f"{len(candidates)} environmental exchange reactions; expected one."
            )
        previous_exchange = environmental_by_nutrient.setdefault(
            nutrient, candidates[0]
        )
        if previous_exchange != candidates[0]:
            raise ValueError(
                f"Nutrient {nutrient!r} maps to multiple ReFramed environmental "
                "exchange reactions."
            )

    reverse = {}
    for nutrient, reaction_id in environmental_by_nutrient.items():
        previous = reverse.setdefault(reaction_id, nutrient)
        if previous != nutrient:
            raise ValueError(
                f"ReFramed exchange {reaction_id!r} collides between nutrients "
                f"{previous!r} and {nutrient!r}."
            )
    if not environmental_by_nutrient:
        raise ValueError(
            "No selected nutrient maps to an actual ReFramed environmental exchange."
        )

    budget_by_reaction = {
        reaction_id: medium[nutrient]
        for nutrient, reaction_id in environmental_by_nutrient.items()
    }
    constraints = {}
    for reaction_id in sorted(exchange_ids):
        reaction = merged.reactions[reaction_id]
        if len(reaction.stoichiometry) != 1:
            raise ValueError(
                f"ReFramed environmental exchange {reaction_id!r} is not a "
                "one-metabolite exchange."
            )
        coefficient = float(next(iter(reaction.stoichiometry.values())))
        if not math.isclose(abs(coefficient), 1.0, rel_tol=0.0, abs_tol=1e-12):
            raise ValueError(
                f"ReFramed environmental exchange {reaction_id!r} has unsupported "
                f"stoichiometric coefficient {coefficient}."
            )
        budget = budget_by_reaction.get(reaction_id, 0.0)
        if coefficient < 0:
            constraints[reaction_id] = (-budget, reaction.ub)
        else:
            constraints[reaction_id] = (reaction.lb, budget)
    return constraints, environmental_by_nutrient


def _validate_solution(solution, model_ids, backend_to_raw):
    growth = float(solution.growth)
    if not math.isfinite(growth) or growth <= 0:
        raise RuntimeError(
            f"SteadyCom did not return finite positive community growth: {growth!r}."
        )
    if set(solution.abundance) != set(backend_to_raw):
        raise RuntimeError(
            "SteadyCom abundance results do not contain the complete model cohort."
        )
    abundance_rows = []
    abundance_sum = 0.0
    zero_ids = []
    for backend, raw_id in backend_to_raw.items():
        abundance = float(solution.abundance[backend])
        if not math.isfinite(abundance) or abundance < 0:
            raise RuntimeError(
                f"SteadyCom returned invalid abundance {abundance!r} for MAG "
                f"{raw_id!r}."
            )
        abundance_sum += abundance
        if abundance == 0:
            zero_ids.append(raw_id)
        abundance_rows.append(
            {"mag": raw_id, "abundance": abundance, "synthetic": False}
        )
    if not math.isclose(abundance_sum, 1.0, rel_tol=1e-6, abs_tol=1e-6):
        raise RuntimeError(
            f"SteadyCom abundances sum to {abundance_sum}, not one within tolerance."
        )
    if zero_ids:
        _log(
            "SteadyCom returned exact zero abundance without pseudocounts for: "
            f"{zero_ids}."
        )

    exchange_rows = []
    for reaction_id, value in solution.exchange.items():
        flux = float(value)
        if not math.isfinite(flux):
            raise RuntimeError(
                f"SteadyCom returned non-finite environmental exchange flux for "
                f"{reaction_id!r}."
            )
        exchange_rows.append(
            {
                "scope": "environment",
                "mag": "",
                "metabolite_id": "",
                "reaction_id": reaction_id,
                "flux": flux,
                "synthetic": False,
            }
        )
    for (backend, raw_species), value in solution.exchange_map.items():
        if backend not in backend_to_raw:
            raise RuntimeError(
                f"SteadyCom exchange_map returned unknown backend ID {backend!r}."
            )
        flux = float(value)
        if not math.isfinite(flux):
            raise RuntimeError(
                f"SteadyCom returned non-finite member exchange rate for "
                f"{backend!r}/{raw_species!r}."
            )
        exchange_rows.append(
            {
                "scope": "member",
                "mag": backend_to_raw[backend],
                "metabolite_id": raw_species,
                "reaction_id": "",
                "flux": flux,
                "synthetic": False,
            }
        )
    if not exchange_rows:
        raise RuntimeError("SteadyCom returned no available exchange values.")
    return growth, abundance_rows, exchange_rows, zero_ids


def _validate_variability_bound(value, raw_id, direction):
    if value is None:
        return None
    try:
        number = float(value)
    except (TypeError, ValueError, OverflowError) as error:
        raise RuntimeError(
            f"SteadyComVA returned an invalid abundance interval: {direction} "
            f"bound {value!r} for MAG {raw_id!r} is not numeric."
        ) from error
    if not math.isfinite(number) or number < -1e-8 or number > 1.0 + 1e-8:
        raise RuntimeError(
            f"SteadyComVA returned an invalid abundance interval: {direction} "
            f"bound {number!r} for MAG {raw_id!r} is outside the finite "
            "abundance range [0, 1]."
        )
    return number


def _validate_variability_interval(lower, upper, raw_id):
    lower = _validate_variability_bound(lower, raw_id, "minimum")
    upper = _validate_variability_bound(upper, raw_id, "maximum")
    if lower is not None and upper is not None and lower > upper + 1e-8:
        raise RuntimeError(
            f"SteadyComVA returned a reversed abundance interval for MAG "
            f"{raw_id!r}: minimum bound {lower} exceeds maximum bound {upper}."
        )
    return lower, upper


def _native_solution_diagnostic(solution):
    status = getattr(solution, "status", None)
    if status is None:
        status_name = None
    elif getattr(status, "name", None) is not None:
        status_name = str(status.name)
    else:
        status_name = str(status)
    message = getattr(solution, "message", None)
    return {
        "native_status": status_name,
        "native_message": None if message is None else str(message),
    }


class _VariabilitySolverObserver:
    """Delegate native solves while retaining SteadyComVA bound diagnostics."""

    def __init__(self, solver, organism_ids, reference_objective=None):
        self._solver = solver
        self._organism_by_variable = {
            f"x_{organism_id}": organism_id for organism_id in organism_ids
        }
        self._reference_objective = reference_objective
        self.reference_diagnostic = None
        self.bound_diagnostics = {}

    def __getattr__(self, name):
        return getattr(self._solver, name)

    def solve(self, linear=None, *args, **kwargs):
        solution = self._solver.solve(linear, *args, **kwargs)
        if isinstance(linear, dict) and len(linear) == 1:
            variable, coefficient = next(iter(linear.items()))
            if (
                variable == self._reference_objective
                and coefficient == 1
                and kwargs.get("get_values", True) is not False
            ):
                self.reference_diagnostic = {
                    **_native_solution_diagnostic(solution),
                    "values_available": getattr(solution, "values", None) is not None,
                }
            organism_id = self._organism_by_variable.get(variable)
            minimize = kwargs.get("minimize")
            if (
                organism_id is not None
                and coefficient == 1
                and isinstance(minimize, bool)
            ):
                direction = "minimum" if minimize else "maximum"
                self.bound_diagnostics[(organism_id, direction)] = (
                    _native_solution_diagnostic(solution)
                )
        return solution


def _missing_bound_diagnostic(raw_id, direction, native_diagnostic):
    if native_diagnostic is None:
        native_diagnostic = {"native_status": None, "native_message": None}
    native_status = native_diagnostic.get("native_status")
    if native_status == "INFEASIBLE":
        bound_status = "infeasible"
    elif native_status == "SUBOPTIMAL":
        bound_status = "nonoptimal"
    elif native_status in (None, "UNKNOWN"):
        bound_status = "unknown"
    else:
        bound_status = "unavailable"
    return {
        "mag": raw_id,
        "direction": direction,
        "bound_status": bound_status,
        "native_status": native_status,
        "native_message": native_diagnostic.get("native_message"),
    }


def _prepare_variability_results(
    variability, backend_to_raw, fraction, bound_diagnostics=None
):
    bound_diagnostics = bound_diagnostics or {}
    rows = []
    missing_bounds = []
    for backend, raw_id in backend_to_raw.items():
        lower, upper = _validate_variability_interval(*variability[backend], raw_id)
        for direction, value in (("minimum", lower), ("maximum", upper)):
            if value is None:
                missing_bounds.append(
                    _missing_bound_diagnostic(
                        raw_id,
                        direction,
                        bound_diagnostics.get((backend, direction)),
                    )
                )
        rows.append(
            {
                "mag": raw_id,
                "minimum_abundance": lower,
                "maximum_abundance": upper,
                "objective_fraction": fraction,
                "variability_status": (
                    "partial" if lower is None or upper is None else "complete"
                ),
                "synthetic": False,
            }
        )
    return rows, missing_bounds, "partial" if missing_bounds else "complete"


def _partial_variability_warning(missing_bounds):
    details = "; ".join(
        f"MAG={item['mag']!r} direction={item['direction']} "
        f"bound_status={item['bound_status']} "
        f"native_status={item['native_status']!r} "
        f"native_message={item['native_message']!r}"
        for item in missing_bounds
    )
    return (
        "WARNING: PARTIAL SteadyComVA results: one or more native variability "
        f"bounds are unavailable ({details}). Missing bounds are exported as NA; "
        "valid SteadyCom growth, abundances, exchanges, and available bounds are "
        "retained."
    )


def _unavailable_variability_results(backend_to_raw, fraction):
    return [
        {
            "mag": raw_id,
            "minimum_abundance": None,
            "maximum_abundance": None,
            "objective_fraction": fraction,
            "variability_status": "unavailable",
            "synthetic": False,
        }
        for raw_id in backend_to_raw.values()
    ]


def _variability_initialization_status(observer, error=None):
    diagnostic = observer.reference_diagnostic or {
        "native_status": None,
        "native_message": None,
        "values_available": False,
    }
    return {
        "result_status": "complete" if error is None else "unavailable",
        **diagnostic,
        "exception_type": None if error is None else type(error).__name__,
        "exception_message": None if error is None else str(error),
    }


def _variability_initialization_warning(initialization_status):
    return (
        "WARNING: PARTIAL SteadyCom results: validated community growth, "
        "abundances, and exchanges are retained, but SteadyComVA initialization "
        "failed and all variability bounds are NOT COMPUTED/unavailable "
        f"(native_status={initialization_status['native_status']!r}, "
        f"native_message={initialization_status['native_message']!r}, "
        f"values_available={initialization_status['values_available']}, "
        f"exception={initialization_status['exception_type']}: "
        f"{initialization_status['exception_message']})."
    )


def _run_variability_analysis(
    steadycom_va,
    community,
    fraction,
    constraints,
    native_solver,
    backend_to_raw,
    log=_log,
):
    reference_objective = community.merged_model.biomass_reaction
    observer = _VariabilitySolverObserver(
        native_solver,
        community.organisms,
        reference_objective=reference_objective,
    )
    try:
        variability = steadycom_va(
            community,
            obj_frac=fraction,
            constraints=constraints,
            solver=observer,
        )
    except Exception as error:
        if observer.bound_diagnostics:
            raise
        initialization_status = _variability_initialization_status(observer, error)
        log(_variability_initialization_warning(initialization_status))
        return (
            _unavailable_variability_results(backend_to_raw, fraction),
            [],
            "unavailable",
            initialization_status,
        )

    if set(variability) != set(backend_to_raw):
        raise RuntimeError(
            "SteadyComVA results do not contain the complete model cohort."
        )
    rows, missing_bounds, variability_status = _prepare_variability_results(
        variability,
        backend_to_raw,
        fraction,
        observer.bound_diagnostics,
    )
    initialization_status = _variability_initialization_status(observer)
    if variability_status == "partial":
        log(_partial_variability_warning(missing_bounds))
    return rows, missing_bounds, variability_status, initialization_status


def _publish_results(
    output,
    growth,
    mock,
    abundance_rows,
    exchange_rows,
    variability_rows,
    status_payload,
):
    _write_tsv(
        output.community_growth,
        [{"stage": "steadycom", "growth": growth, "synthetic": mock}],
        ["stage", "growth", "synthetic"],
    )
    _write_tsv(
        output.abundances,
        abundance_rows,
        ["mag", "abundance", "synthetic"],
    )
    _write_tsv(
        output.exchanges,
        exchange_rows,
        ["scope", "mag", "metabolite_id", "reaction_id", "flux", "synthetic"],
    )
    _write_tsv(
        output.variability,
        variability_rows,
        [
            "mag",
            "minimum_abundance",
            "maximum_abundance",
            "objective_fraction",
            "variability_status",
            "synthetic",
        ],
        na_rep="NA",
    )
    Path(output.status).parent.mkdir(parents=True, exist_ok=True)
    Path(output.status).write_text(
        json.dumps(status_payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def main():
    tag = str(snakemake.params.tag)
    solver_name = str(snakemake.params.solver)
    objective = str(snakemake.params.objective)
    fraction = float(snakemake.params.variability_fraction)
    mock = bool(snakemake.params.mock)
    model_ids, model_paths = model_ids_and_paths(
        snakemake.input.models, snakemake.input.filtered_models
    )
    medium = read_medium_csv(snakemake.input.medium)
    medium_status = _validate_status(snakemake.input.medium_status, tag)
    raw_to_backend = backend_ids(model_ids)
    backend_to_raw = {backend: raw for raw, backend in raw_to_backend.items()}
    if objective != "community_growth":
        raise ValueError(
            "The pinned native SteadyCom path supports objective=community_growth only."
        )

    if mock:
        growth = 0.0
        abundance_rows = [
            {
                "mag": model_id,
                "abundance": 1.0 / len(model_ids),
                "synthetic": True,
            }
            for model_id in model_ids
        ]
        exchange_rows = [
            {
                "scope": "environment",
                "mag": "",
                "metabolite_id": nutrient,
                "reaction_id": f"synthetic_import_{nutrient}",
                "flux": -uptake,
                "synthetic": True,
            }
            for nutrient, uptake in medium.items()
        ]
        variability_rows = [
            {
                "mag": model_id,
                "minimum_abundance": 0.0,
                "maximum_abundance": 1.0,
                "objective_fraction": fraction,
                "variability_status": "complete",
                "synthetic": True,
            }
            for model_id in model_ids
        ]
        environmental_map = {
            nutrient: f"synthetic_import_{nutrient}" for nutrient in medium
        }
        zero_ids = []
        missing_variability_bounds = []
        variability_status = "complete"
        variability_initialization = {
            "result_status": "complete",
            "native_status": None,
            "native_message": None,
            "values_available": True,
            "exception_type": None,
            "exception_message": None,
        }
    else:
        if solver_name != "cplex":
            raise ValueError(
                "Native ReFramed SteadyCom/SteadyComVA requires solver=cplex; "
                "the pinned growth-update implementation is CPLEX-specific."
            )
        ensure_cplex_ready(
            solver_name,
            snakemake.config.get("cplex", {}).get("studio_dir", ""),
            auto_upgrade=bool(
                snakemake.config.get("cplex", {}).get("auto_upgrade", True)
            ),
            log_file=str(snakemake.log[0]),
        )
        try:
            from reframed import Community, load_cbmodel, set_default_solver
            from reframed.community.SteadyCom import (
                SteadyCom,
                SteadyComVA,
                build_problem,
            )
        except Exception as error:
            raise RuntimeError(
                "SteadyCom requires ReFramed 1.5.1 with native SteadyCom and "
                "SteadyComVA support."
            ) from error

        set_default_solver(solver_name)
        organisms = []
        for model_id, model_path in zip(model_ids, model_paths):
            organism = load_cbmodel(model_path, flavor="fbc2")
            organism.id = raw_to_backend[model_id]
            organisms.append(organism)
        community = Community("community", organisms)
        if set(community.organisms) != set(backend_to_raw):
            raise RuntimeError(
                "ReFramed did not retain the complete model cohort; no models "
                "may be silently dropped."
            )
        constraints, environmental_map = _medium_constraints(
            community,
            read_mapping_tsv(snakemake.input.medium_mapping),
            medium,
            raw_to_backend,
        )
        native_solver = build_problem(community)
        native_solver.problem.parameters.threads.set(int(snakemake.threads))
        _log(
            f"Running ReFramed SteadyCom and SteadyComVA for {len(model_ids)} "
            f"model(s) on selected medium {tag!r}."
        )
        solution = SteadyCom(community, constraints=constraints, solver=native_solver)
        growth, abundance_rows, exchange_rows, zero_ids = _validate_solution(
            solution, model_ids, backend_to_raw
        )
        (
            variability_rows,
            missing_variability_bounds,
            variability_status,
            variability_initialization,
        ) = _run_variability_analysis(
            SteadyComVA,
            community,
            fraction,
            constraints,
            native_solver,
            backend_to_raw,
        )

    analysis_status = "complete" if variability_status == "complete" else "partial"
    status_payload = {
        "synthetic": mock,
        "tool": "steadycom",
        "medium_tag": tag,
        "medium_source": medium_status.get("source"),
        "solver": solver_name,
        "threads": int(snakemake.threads),
        "objective": objective,
        "variability_fraction": fraction,
        "analysis_status": analysis_status,
        "variability_status": variability_status,
        "stage_status": {
            "steadycom": {"result_status": "complete"},
            "steadycom_va": {
                "result_status": variability_status,
                "initialization": variability_initialization,
            },
        },
        "missing_variability_bounds": missing_variability_bounds,
        "model_ids": model_ids,
        "zero_abundance_model_ids": zero_ids,
        "environmental_exchange_by_nutrient": environmental_map,
    }
    _publish_results(
        snakemake.output,
        growth,
        mock,
        abundance_rows,
        exchange_rows,
        variability_rows,
        status_payload,
    )


if "snakemake" in globals():
    try:
        main()
    except Exception:
        _log("SteadyCom failed:\n" + traceback.format_exc())
        raise
