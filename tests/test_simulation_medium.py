import copy
import csv
import json
import math
import sys
from enum import Enum
from pathlib import Path
from types import SimpleNamespace

import pytest
import yaml
from jsonschema import Draft202012Validator, ValidationError

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

from common.media import (  # noqa: E402
    aggregate_media,
    backend_ids,
    model_ids_and_paths,
    read_medium_csv,
)
from micom_script import _validate_solution_metrics  # noqa: E402
from steadycom import (  # noqa: E402
    _partial_variability_warning,
    _prepare_variability_results,
    _publish_results,
    _run_variability_analysis,
    _validate_variability_interval,
    _VariabilitySolverObserver,
)


def _medium(tmp_path, name, rows):
    path = tmp_path / name
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["nutrient_id", "uptake"])
        writer.writerows(rows)
    return path


def test_sum_default_and_max_include_missing_nutrients_as_zero(tmp_path):
    first = read_medium_csv(_medium(tmp_path, "first.csv", [("a_e", 2), ("b_e", 3)]))
    second = read_medium_csv(_medium(tmp_path, "second.csv", [("a_e", 5), ("c_e", 7)]))

    assert aggregate_media([first, second], "sum") == {
        "a_e": 7.0,
        "b_e": 3.0,
        "c_e": 7.0,
    }
    assert aggregate_media([first, second], "max") == {
        "a_e": 5.0,
        "b_e": 3.0,
        "c_e": 7.0,
    }


@pytest.mark.parametrize(
    ("rows", "message"),
    [
        ([("a_e", 1), ("a_e", 2)], "duplicate"),
        ([("a_e", 0)], "finite and positive"),
        ([("a_e", -1)], "finite and positive"),
        ([("a_e", math.inf)], "finite and positive"),
        ([("", 1)], "empty nutrient_id"),
        ([("a_e", 1, "extra")], "exactly two values"),
    ],
)
def test_medium_csv_rejects_duplicate_nonfinite_and_nonpositive_values(
    tmp_path, rows, message
):
    with pytest.raises(ValueError, match=message):
        read_medium_csv(_medium(tmp_path, "bad.csv", rows))


def test_m9_is_exact_approved_20_by_100_with_oxygen():
    medium = read_medium_csv(ROOT / "workflow" / "resources" / "M9.csv")
    assert len(medium) == 20
    assert set(medium) == {
        "ca2_e",
        "cl_e",
        "cobalt2_e",
        "cu2_e",
        "fe2_e",
        "fe3_e",
        "glc__D_e",
        "h2o_e",
        "h_e",
        "k_e",
        "mg2_e",
        "mn2_e",
        "mobd_e",
        "na1_e",
        "nh4_e",
        "ni2_e",
        "o2_e",
        "pi_e",
        "so4_e",
        "zn2_e",
    }
    assert set(medium.values()) == {100.0}


def test_complete_retained_cohort_and_reversible_backend_identity(tmp_path):
    models = [tmp_path / "raw-A.xml", tmp_path / "raw.B.xml"]
    for model in models:
        model.touch()
    retained = tmp_path / "filtered_models.txt"
    retained.write_text("raw.B\nraw-A\n", encoding="utf-8")

    ids, ordered = model_ids_and_paths(models, retained)
    assert ids == ["raw.B", "raw-A"]
    assert list(map(Path, ordered)) == [models[1], models[0]]
    mapping = backend_ids(ids)
    assert mapping == {"raw.B": "member_0001", "raw-A": "member_0002"}
    assert {backend: raw for raw, backend in mapping.items()} == {
        "member_0001": "raw.B",
        "member_0002": "raw-A",
    }


def test_missing_or_duplicate_retained_recipe_identity_is_an_error(tmp_path):
    model = tmp_path / "one.xml"
    model.touch()
    retained = tmp_path / "filtered_models.txt"
    retained.write_text("one\ntwo\n", encoding="utf-8")
    with pytest.raises(ValueError, match="missing model files"):
        model_ids_and_paths([model], retained)

    retained.write_text("one\none\n", encoding="utf-8")
    with pytest.raises(ValueError, match="duplicate MAG"):
        model_ids_and_paths([model], retained)


def test_simulation_medium_schema_contract_and_defaults():
    schema = yaml.safe_load(
        (ROOT / "workflow" / "schemas" / "config.schema.yaml").read_text()
    )
    config = yaml.safe_load((ROOT / "config" / "config.yaml").read_text())
    validator = Draft202012Validator(schema)
    validator.validate(config)
    assert config["simulation_medium"] == {
        "source": "discover",
        "aggregation": "sum",
        "predefined_medium": "M9",
        "min_growth": 0.1,
        "open_exchanges": True,
        "exchange_bound": 1000,
    }

    for key, value in (
        ("source", "community"),
        ("aggregation", "mean"),
        ("min_growth", 0),
        ("exchange_bound", "unbounded"),
    ):
        invalid = copy.deepcopy(config)
        invalid["simulation_medium"][key] = value
        with pytest.raises(ValidationError):
            validator.validate(invalid)


@pytest.mark.parametrize("status", ["suboptimal", "optimal_inaccurate", "OPTIMAL"])
def test_micom_rejects_nonexact_optimal_status(status):
    with pytest.raises(RuntimeError, match="finite optimal solution"):
        _validate_solution_metrics("fba", status, 0.1, 0.1)


def test_steadycom_variability_rejects_upper_bound_above_one():
    with pytest.raises(RuntimeError, match="invalid abundance interval"):
        _validate_variability_interval(0.2, 1.00000002, "mag-a")

    assert _validate_variability_interval(0.0, 1.000000005, "mag-a") == (
        0.0,
        1.000000005,
    )


@pytest.mark.parametrize(
    ("lower", "upper", "missing_directions"),
    [
        (None, 0.8, ["minimum"]),
        (0.2, None, ["maximum"]),
        (None, None, ["minimum", "maximum"]),
    ],
    ids=["minimum-none", "maximum-none", "both-none"],
)
def test_steadycom_variability_none_bounds_are_explicit_partial_results(
    lower, upper, missing_directions
):
    diagnostics = {
        ("member_0001", "minimum"): {
            "native_status": "INFEASIBLE",
            "native_message": "infeasible",
        },
        ("member_0001", "maximum"): {
            "native_status": "UNKNOWN",
            "native_message": "numerical difficulty",
        },
    }

    rows, missing, overall = _prepare_variability_results(
        {"member_0001": (lower, upper)},
        {"member_0001": "mag-a"},
        0.9,
        diagnostics,
    )

    assert overall == "partial"
    assert rows == [
        {
            "mag": "mag-a",
            "minimum_abundance": lower,
            "maximum_abundance": upper,
            "objective_fraction": 0.9,
            "variability_status": "partial",
            "synthetic": False,
        }
    ]
    assert [item["direction"] for item in missing] == missing_directions
    for item in missing:
        assert item["mag"] == "mag-a"
        if item["direction"] == "minimum":
            assert item == {
                "mag": "mag-a",
                "direction": "minimum",
                "bound_status": "infeasible",
                "native_status": "INFEASIBLE",
                "native_message": "infeasible",
            }
        else:
            assert item == {
                "mag": "mag-a",
                "direction": "maximum",
                "bound_status": "unknown",
                "native_status": "UNKNOWN",
                "native_message": "numerical difficulty",
            }


def test_steadycom_variability_preserves_genuine_zero_and_complete_results():
    rows, missing, overall = _prepare_variability_results(
        {"member_0001": (0.0, 0.75)},
        {"member_0001": "mag-zero"},
        0.9,
    )

    assert overall == "complete"
    assert missing == []
    assert rows[0]["minimum_abundance"] == 0.0
    assert rows[0]["maximum_abundance"] == 0.75
    assert rows[0]["variability_status"] == "complete"


@pytest.mark.parametrize(
    ("lower", "upper", "directions"),
    [
        (math.nan, 0.5, ("minimum",)),
        (0.1, math.inf, ("maximum",)),
        (-0.001, 0.5, ("minimum",)),
        (0.1, 1.001, ("maximum",)),
        ("not-a-number", 0.5, ("minimum",)),
        (0.8, 0.2, ("minimum", "maximum")),
        (None, 1.001, ("maximum",)),
    ],
    ids=[
        "nan",
        "infinite",
        "lower-out-of-range",
        "upper-out-of-range",
        "nonnumeric",
        "reversed",
        "missing-with-invalid-present-bound",
    ],
)
def test_steadycom_variability_invalid_present_bounds_still_fail_with_context(
    lower, upper, directions
):
    with pytest.raises(RuntimeError) as caught:
        _validate_variability_interval(lower, upper, "mag-invalid")

    message = str(caught.value)
    assert "mag-invalid" in message
    assert all(direction in message for direction in directions)


class _InjectedStatus(Enum):
    OPTIMAL = "Optimal"
    INFEASIBLE = "Infeasible"
    UNKNOWN = "Unknown"
    SUBOPTIMAL = "Suboptimal"


class _InjectedSolution:
    def __init__(self, status, message, fobj=None, values=None):
        self.status = status
        self.message = message
        self.fobj = fobj
        self.values = values


class _InjectedSolver:
    def __init__(self, solutions):
        self.solutions = iter(solutions)
        self.calls = []

    def solve(self, linear=None, *args, **kwargs):
        self.calls.append((linear, args, kwargs))
        return next(self.solutions)


def test_injected_solver_statuses_are_observed_without_none_inference():
    solver = _InjectedSolver(
        [
            _InjectedSolution(_InjectedStatus.INFEASIBLE, "infeasible"),
            _InjectedSolution(_InjectedStatus.SUBOPTIMAL, "stopped early"),
        ]
    )
    observer = _VariabilitySolverObserver(solver, ["member_0001"])

    minimum = observer.solve({"x_member_0001": 1}, minimize=True, get_values=False)
    maximum = observer.solve({"x_member_0001": 1}, minimize=False, get_values=False)

    assert minimum.fobj is None and maximum.fobj is None
    assert observer.bound_diagnostics == {
        ("member_0001", "minimum"): {
            "native_status": "INFEASIBLE",
            "native_message": "infeasible",
        },
        ("member_0001", "maximum"): {
            "native_status": "SUBOPTIMAL",
            "native_message": "stopped early",
        },
    }
    _, missing, _ = _prepare_variability_results(
        {"member_0001": (None, None)},
        {"member_0001": "mag-a"},
        0.9,
        observer.bound_diagnostics,
    )
    assert [item["bound_status"] for item in missing] == [
        "infeasible",
        "nonoptimal",
    ]
    assert solver.calls == [
        ({"x_member_0001": 1}, (), {"minimize": True, "get_values": False}),
        ({"x_member_0001": 1}, (), {"minimize": False, "get_values": False}),
    ]


def _injected_community():
    return SimpleNamespace(
        organisms={"member_0001": object(), "member_0002": object()},
        merged_model=SimpleNamespace(biomass_reaction="community_biomass"),
    )


def test_injected_va_missing_reference_values_is_unavailable_not_bound_infeasible():
    solver = _InjectedSolver(
        [
            _InjectedSolution(
                _InjectedStatus.INFEASIBLE,
                "infeasible",
                values=None,
            )
        ]
    )
    warnings = []

    def missing_reference_values(community, obj_frac, constraints, solver):
        reference = solver.solve({"community_biomass": 1})
        return reference.values["community_biomass"]

    rows, missing, status, initialization = _run_variability_analysis(
        missing_reference_values,
        _injected_community(),
        0.9,
        {},
        solver,
        {"member_0001": "mag-a", "member_0002": "mag-b"},
        log=warnings.append,
    )

    assert status == "unavailable"
    assert missing == []
    assert all(row["minimum_abundance"] is None for row in rows)
    assert all(row["maximum_abundance"] is None for row in rows)
    assert all(row["variability_status"] == "unavailable" for row in rows)
    assert initialization == {
        "result_status": "unavailable",
        "native_status": "INFEASIBLE",
        "native_message": "infeasible",
        "values_available": False,
        "exception_type": "TypeError",
        "exception_message": "'NoneType' object is not subscriptable",
    }
    assert "NOT COMPUTED/unavailable" in warnings[0]
    assert "growth, abundances, and exchanges are retained" in warnings[0]


def test_injected_va_optimizer_exception_is_unavailable_with_unknown_status():
    warnings = []

    def optimizer_exception(community, obj_frac, constraints, solver):
        raise RuntimeError("reference solve exploded")

    rows, missing, status, initialization = _run_variability_analysis(
        optimizer_exception,
        _injected_community(),
        0.9,
        {},
        _InjectedSolver([]),
        {"member_0001": "mag-a", "member_0002": "mag-b"},
        log=warnings.append,
    )

    assert status == "unavailable" and missing == []
    assert len(rows) == 2
    assert initialization["native_status"] is None
    assert initialization["exception_type"] == "RuntimeError"
    assert initialization["exception_message"] == "reference solve exploded"


def test_injected_va_ignores_normal_nonoptimal_trial_and_records_final_reference():
    solver = _InjectedSolver(
        [
            _InjectedSolution(_InjectedStatus.INFEASIBLE, "trial infeasible"),
            _InjectedSolution(
                _InjectedStatus.OPTIMAL,
                "optimal",
                values={"community_biomass": 1.0},
            ),
            _InjectedSolution(_InjectedStatus.OPTIMAL, "optimal", fobj=0.0),
            _InjectedSolution(_InjectedStatus.OPTIMAL, "optimal", fobj=0.8),
            _InjectedSolution(_InjectedStatus.OPTIMAL, "optimal", fobj=0.2),
            _InjectedSolution(_InjectedStatus.OPTIMAL, "optimal", fobj=1.0),
        ]
    )

    def successful_va(community, obj_frac, constraints, solver):
        solver.solve({"community_biomass": 1}, get_values=False)
        reference = solver.solve({"community_biomass": 1})
        assert reference.values["community_biomass"] == 1.0
        return {
            "member_0001": [
                solver.solve(
                    {"x_member_0001": 1}, minimize=True, get_values=False
                ).fobj,
                solver.solve(
                    {"x_member_0001": 1}, minimize=False, get_values=False
                ).fobj,
            ],
            "member_0002": [
                solver.solve(
                    {"x_member_0002": 1}, minimize=True, get_values=False
                ).fobj,
                solver.solve(
                    {"x_member_0002": 1}, minimize=False, get_values=False
                ).fobj,
            ],
        }

    rows, missing, status, initialization = _run_variability_analysis(
        successful_va,
        _injected_community(),
        0.9,
        {},
        solver,
        {"member_0001": "mag-a", "member_0002": "mag-b"},
        log=lambda message: None,
    )

    assert status == "complete" and missing == []
    assert [(row["minimum_abundance"], row["maximum_abundance"]) for row in rows] == [
        (0.0, 0.8),
        (0.2, 1.0),
    ]
    assert initialization["result_status"] == "complete"
    assert initialization["native_status"] == "OPTIMAL"
    assert initialization["values_available"] is True


def test_partial_export_writes_na_warning_status_and_preserves_main_outputs(tmp_path):
    diagnostics = {
        ("member_0001", "maximum"): {
            "native_status": "UNKNOWN",
            "native_message": "numerical difficulty",
        }
    }
    rows, missing, overall = _prepare_variability_results(
        {
            "member_0001": (0.0, None),
            "member_0002": (0.25, 0.75),
        },
        {"member_0001": "mag-a", "member_0002": "mag-b"},
        0.9,
        diagnostics,
    )
    output = SimpleNamespace(
        community_growth=tmp_path / "community_growth.tsv",
        abundances=tmp_path / "abundances.tsv",
        exchanges=tmp_path / "exchanges.tsv",
        variability=tmp_path / "variability.tsv",
        status=tmp_path / "status.json",
    )
    status_payload = {
        "tool": "steadycom",
        "variability_status": overall,
        "missing_variability_bounds": missing,
    }

    _publish_results(
        output,
        0.5,
        False,
        [
            {"mag": "mag-a", "abundance": 0.4, "synthetic": False},
            {"mag": "mag-b", "abundance": 0.6, "synthetic": False},
        ],
        [
            {
                "scope": "environment",
                "mag": "",
                "metabolite_id": "",
                "reaction_id": "R_EX_glc",
                "flux": -0.5,
                "synthetic": False,
            }
        ],
        rows,
        status_payload,
    )

    assert all(
        path.is_file()
        for path in (
            output.community_growth,
            output.abundances,
            output.exchanges,
            output.variability,
            output.status,
        )
    )
    with output.variability.open(newline="", encoding="utf-8") as handle:
        exported = list(csv.DictReader(handle, delimiter="\t"))
    assert exported[0]["minimum_abundance"] == "0.0"
    assert exported[0]["maximum_abundance"] == "NA"
    assert exported[0]["variability_status"] == "partial"
    assert exported[1]["minimum_abundance"] == "0.25"
    assert exported[1]["maximum_abundance"] == "0.75"
    assert exported[1]["variability_status"] == "complete"
    assert "\t0.5\tFalse\n" in output.community_growth.read_text(encoding="utf-8")
    assert "mag-a\t0.4\tFalse" in output.abundances.read_text(encoding="utf-8")
    assert "R_EX_glc\t-0.5\tFalse" in output.exchanges.read_text(encoding="utf-8")
    assert json.loads(output.status.read_text(encoding="utf-8")) == status_payload
    warning = _partial_variability_warning(missing)
    assert "WARNING: PARTIAL SteadyComVA results" in warning
    assert "Missing bounds are exported as NA" in warning
    assert "native_status='UNKNOWN'" in warning
