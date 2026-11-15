"""Opt-in native regression for the selected simulation-medium workflow."""

import csv
import json
import os
import subprocess
import time
import uuid
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
SNAKEFILE = ROOT / "workflow" / "Snakefile"
FIXTURE = ROOT / "tests" / "native_medium_fixture.py"
MAIN_ROOT = Path(os.environ.get("MEEW_MAIN_ROOT", ROOT.parent / "MeeW"))
MAIN_SNAKEMAKE = MAIN_ROOT / ".snakemake"

pytestmark = pytest.mark.skipif(
    os.environ.get("MEEW_RUN_NATIVE_MEDIUM_TESTS") != "1",
    reason="set MEEW_RUN_NATIVE_MEDIUM_TESTS=1 for pinned CPLEX native tests",
)


def _rows(path, delimiter=","):
    with Path(path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def _medium(path):
    return {row["nutrient_id"]: float(row["uptake"]) for row in _rows(path)}


def _write_config(root, **medium):
    config = yaml.safe_load((ROOT / "config" / "config.yaml").read_text())
    config["input"]["faa_globs"] = ["data/*.faa"]
    config["sample"] = "native"
    config["qa"]["enabled"] = False
    config["checkm2"]["enabled"] = False
    config["memote"]["min_total_score"] = 0.0
    config["community"] = {
        "enabled": True,
        "tools": {"smetana": False, "micom": True, "steadycom": True},
    }
    config["simulation_medium"].update(medium)
    config["micom"]["abundance_tsv"] = (
        str(root / "data" / "abundance.tsv")
        if (root / "scenario.txt").read_text(encoding="utf-8") == "main"
        else ""
    )
    config["execution"]["mock"] = False
    path = root / "config.yaml"
    path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    return path


def _change_medium(config_path, **updates):
    time.sleep(0.02)
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    config["simulation_medium"].update(updates)
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")


def _new_case(parent, name, scenario, micom_python):
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    root = parent / f"{stamp}-{name}-{uuid.uuid4().hex[:8]}"
    root.mkdir(parents=True)
    (root / "scenario.txt").write_text(scenario, encoding="utf-8")
    generated = subprocess.run(
        [
            str(micom_python),
            str(FIXTURE),
            "--root",
            str(root),
            "--scenario",
            scenario,
        ],
        cwd=ROOT,
        env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
        check=False,
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert generated.returncode == 0, generated.stdout + generated.stderr
    return root


def _run(root, config_path, conda_prefix, label, *targets):
    command = [
        "snakemake",
        "--snakefile",
        str(SNAKEFILE),
        "--directory",
        str(root),
        "--configfile",
        str(config_path),
        "--replace-workflow-config",
        "--nolock",
        "--cores",
        "2",
        "--rerun-triggers=mtime",
        "--use-conda",
        "--conda-prefix",
        str(conda_prefix),
        "--printshellcmds",
        "--show-failed-logs",
        *targets,
    ]
    result = subprocess.run(
        command,
        cwd=ROOT,
        env={**os.environ, "GIT_OPTIONAL_LOCKS": "0", "PYTHONDONTWRITEBYTECODE": "1"},
        check=False,
        capture_output=True,
        text=True,
        timeout=240,
    )
    logs = root / "run-logs"
    logs.mkdir(exist_ok=True)
    (logs / f"{label}.command.txt").write_text(
        " ".join(command) + "\n", encoding="utf-8"
    )
    (logs / f"{label}.stdout.txt").write_text(result.stdout, encoding="utf-8")
    (logs / f"{label}.stderr.txt").write_text(result.stderr, encoding="utf-8")
    return result


def _assert_success(result):
    assert result.returncode == 0, result.stdout + result.stderr
    return result.stdout + result.stderr


def _source_snapshot(root):
    return {
        path: (path.read_bytes(), path.stat().st_mtime_ns)
        for path in sorted((root / "results" / "native" / "gems").glob("*.xml"))
    }


def _assert_source_unchanged(snapshot):
    assert {
        path: (path.read_bytes(), path.stat().st_mtime_ns) for path in snapshot
    } == snapshot


def test_native_cplex_genuine_infeasible_status_is_observed():
    steady_python = Path(os.environ["MEEW_NATIVE_STEADYCOM_PYTHON"]).resolve()
    script = f"""
import json
import sys
sys.path.insert(0, {str(ROOT / 'workflow' / 'scripts')!r})
from reframed import set_default_solver, solver_instance
from steadycom import _VariabilitySolverObserver

set_default_solver('cplex')
solver = solver_instance()
solver.add_variable('x_member_0001', 0, 1, update=False)
solver.add_constraint(
    'genuinely_impossible', {{'x_member_0001': 1}}, '>', 2, update=False
)
solver.update()
observer = _VariabilitySolverObserver(solver, ['member_0001'])
solution = observer.solve(
    {{'x_member_0001': 1}}, minimize=True, get_values=False
)
print(json.dumps({{
    'fobj': solution.fobj,
    'diagnostic': observer.bound_diagnostics[
        ('member_0001', 'minimum')
    ],
}}))
"""
    result = subprocess.run(
        [str(steady_python), "-c", script],
        cwd=ROOT,
        env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
        check=False,
        capture_output=True,
        text=True,
        timeout=60,
    )

    assert result.returncode == 0, result.stdout + result.stderr
    payload = json.loads([line for line in result.stdout.splitlines() if line][-1])
    assert payload["fobj"] is None
    assert payload["diagnostic"]["native_status"] == "INFEASIBLE"
    assert "infeasible" in payload["diagnostic"]["native_message"].lower()


def test_native_selected_medium_regression():
    native_parent = Path(
        os.environ.get(
            "MEEW_NATIVE_TEST_ROOT", MAIN_SNAKEMAKE / "native-medium-regressions"
        )
    ).resolve()
    main_snakemake = MAIN_SNAKEMAKE.resolve()
    assert native_parent.is_relative_to(
        main_snakemake
    ), f"MEEW_NATIVE_TEST_ROOT must remain under {main_snakemake}"
    conda_prefix = Path(os.environ["MEEW_NATIVE_CONDA_PREFIX"]).resolve()
    micom_python = Path(os.environ["MEEW_NATIVE_MICOM_PYTHON"]).resolve()
    steady_python = Path(os.environ["MEEW_NATIVE_STEADYCOM_PYTHON"]).resolve()
    assert micom_python.is_file() and steady_python.is_file()
    assert micom_python.parents[1].parent == conda_prefix
    assert steady_python.parents[1].parent == conda_prefix
    native_parent.mkdir(parents=True, exist_ok=True)

    main = _new_case(native_parent, "all-modes", "main", micom_python)
    config = _write_config(
        main,
        source="discover",
        aggregation="sum",
        predefined_medium="M9",
        min_growth=0.1,
        open_exchanges=True,
        exchange_bound=1000,
    )
    source_before = _source_snapshot(main)
    results = main / "results" / "native"

    rendered = _assert_success(_run(main, config, conda_prefix, "01-sum"))
    assert "localrule discover_individual_medium:" in rendered
    assert micom_python.parents[1].name in rendered
    assert steady_python.parents[1].name in rendered
    assert _medium(results / "media" / "summed_demand.csv") == pytest.approx(
        {"glc__D_e": 0.3, "o2_e": 0.2}
    )
    mappings = _rows(results / "media" / "summed_demand.mapping.tsv", "\t")
    assert {(row["nutrient_id"], row["cobra_exchange_id"]) for row in mappings} >= {
        ("glc__D_e", "EX_glc__D_e"),
        ("o2_e", "EX_o2_e"),
    }
    source_root = ET.parse(results / "gems" / "raw-A.xml").getroot()
    source_reactions = {
        element.attrib.get("id")
        for element in source_root.iter()
        if element.tag.rsplit("}", 1)[-1] == "reaction"
    }
    source_oxygen = [
        element
        for element in source_root.iter()
        if element.tag.rsplit("}", 1)[-1] == "species"
        and element.attrib.get("id") == "M_o2_e"
    ]
    assert "R_EX_glc__D_e" in source_reactions
    assert "R_EX_o2_e" not in source_reactions
    assert len(source_oxygen) == 1
    assert source_oxygen[0].attrib["boundaryCondition"] == "true"
    for model_id in ("raw-A", "raw.B"):
        summary = json.loads((results / "solo" / model_id / "summary.json").read_text())
        assert summary["verification_status"] == "optimal"
        assert summary["verified_growth"] >= 0.1
        assert summary["threads"] == 1

    micom_status = json.loads(
        (results / "micom" / "summed_demand" / "status.json").read_text()
    )
    assert micom_status["model_ids"] == ["raw-A", "raw.B"]
    assert micom_status["zero_abundance_model_ids"] == ["raw.B"]
    assert micom_status["threads"] == 1
    assert micom_status["analysis_status"] == "complete"
    assert micom_status["stage_status"]["fba"]["result_status"] == "complete"
    assert micom_status["stage_status"]["pfba"]["result_status"] == "complete"
    members = _rows(results / "micom" / "summed_demand" / "member_growth.tsv", "\t")
    assert {(row["stage"], row["mag"]) for row in members} == {
        (stage, mag) for stage in ("fba", "pfba") for mag in ("raw-A", "raw.B")
    }
    assert all(
        float(row["abundance"]) == 0.0 for row in members if row["mag"] == "raw.B"
    )
    assert {row["stage"] for row in members} == {"fba", "pfba"}
    assert any(
        float(row["flux"]) != 0.0
        for row in _rows(results / "micom" / "summed_demand" / "fluxes.tsv", "\t")
    )

    steady_status = json.loads(
        (results / "steadycom" / "summed_demand" / "status.json").read_text()
    )
    assert steady_status["model_ids"] == ["raw-A", "raw.B"]
    assert steady_status["zero_abundance_model_ids"] == ["raw.B"]
    assert steady_status["threads"] == 1
    assert steady_status["analysis_status"] == "complete"
    assert steady_status["variability_status"] == "complete"
    assert steady_status["missing_variability_bounds"] == []
    assert (
        steady_status["stage_status"]["steadycom_va"]["initialization"]["result_status"]
        == "complete"
    )
    abundances = _rows(results / "steadycom" / "summed_demand" / "abundances.tsv", "\t")
    assert {row["mag"]: float(row["abundance"]) for row in abundances} == {
        "raw-A": 1.0,
        "raw.B": 0.0,
    }
    variability = _rows(
        results / "steadycom" / "summed_demand" / "variability.tsv", "\t"
    )
    assert all(float(row["objective_fraction"]) == 0.9 for row in variability)
    assert all(row["variability_status"] == "complete" for row in variability)
    assert all(
        row[bound] != "NA"
        for row in variability
        for bound in ("minimum_abundance", "maximum_abundance")
    )

    cached = _assert_success(_run(main, config, conda_prefix, "02-sum-cached"))
    assert "Nothing to be done" in cached

    _change_medium(config, aggregation="max")
    _assert_success(_run(main, config, conda_prefix, "03-max"))
    assert _medium(results / "media" / "max_demand.csv") == pytest.approx(
        {"glc__D_e": 0.2, "o2_e": 0.1}
    )

    _change_medium(config, source="predefined", predefined_medium="M9")
    m9_run = _assert_success(_run(main, config, conda_prefix, "04-m9"))
    assert "localrule discover_individual_medium:" not in m9_run
    m9 = _medium(results / "media" / "M9.csv")
    assert len(m9) == 20 and set(m9.values()) == {100.0} and m9["o2_e"] == 100

    custom_a = main / "data" / "custom-a.csv"
    custom_b = main / "data" / "custom-b.csv"
    custom_a.write_text("nutrient_id,uptake\nglc__D_e,0.5\no2_e,0.5\n")
    custom_b.write_text("nutrient_id,uptake\nglc__D_e,0.4\no2_e,0.4\n")
    os.utime(custom_b, ns=(1_000_000_000, 1_000_000_000))
    _change_medium(config, predefined_medium=str(custom_a))
    custom_run = _assert_success(_run(main, config, conda_prefix, "05-custom-a"))
    assert "localrule discover_individual_medium:" not in custom_run
    assert _medium(results / "media" / "custom.csv")["glc__D_e"] == 0.5
    assert "Nothing to be done" in _assert_success(
        _run(main, config, conda_prefix, "06-custom-a-cached")
    )
    custom_output_mtime = (results / "media" / "custom.csv").stat().st_mtime_ns
    assert custom_b.stat().st_mtime_ns < custom_output_mtime
    _change_medium(config, predefined_medium=str(custom_b))
    replaced = _assert_success(_run(main, config, conda_prefix, "07-custom-b-older"))
    assert "localrule prepare_simulation_medium:" in replaced
    assert "localrule discover_individual_medium:" not in replaced
    assert _medium(results / "media" / "custom.csv")["glc__D_e"] == 0.4
    for tag in ("summed_demand", "max_demand", "M9", "custom"):
        for tool in ("micom", "steadycom"):
            status = json.loads((results / tool / tag / "status.json").read_text())
            assert status["synthetic"] is False
            assert status["model_ids"] == ["raw-A", "raw.B"]
    _assert_source_unchanged(source_before)

    blocked_open = _new_case(native_parent, "blocked-open", "blocked", micom_python)
    open_config = _write_config(
        blocked_open,
        source="discover",
        aggregation="sum",
        predefined_medium="M9",
        min_growth=0.1,
        open_exchanges=True,
        exchange_bound=1000,
    )
    blocked_source = _source_snapshot(blocked_open)
    open_target = "results/native/solo/blocked/minimal_medium.csv"
    _assert_success(
        _run(blocked_open, open_config, conda_prefix, "01-open", open_target)
    )
    open_summary = json.loads(
        (
            blocked_open / "results" / "native" / "solo" / "blocked" / "summary.json"
        ).read_text()
    )
    assert open_summary["open_exchanges"] is True
    assert open_summary["verification_status"] == "optimal"
    assert open_summary["verified_growth"] >= 0.1
    _assert_source_unchanged(blocked_source)

    blocked_closed = _new_case(native_parent, "blocked-closed", "blocked", micom_python)
    closed_config = _write_config(
        blocked_closed,
        source="discover",
        aggregation="sum",
        predefined_medium="M9",
        min_growth=0.1,
        open_exchanges=False,
        exchange_bound=1000,
    )
    closed_source = _source_snapshot(blocked_closed)
    closed = _run(blocked_closed, closed_config, conda_prefix, "01-closed", open_target)
    assert closed.returncode != 0
    assert "No minimum-import medium supports growth 0.1" in (
        closed.stdout + closed.stderr
    )
    assert (
        "open_exchanges=False"
        in (
            blocked_closed / "results" / "native" / "logs" / "solo" / "blocked.log"
        ).read_text()
    )
    _assert_source_unchanged(closed_source)

    orphan = _new_case(native_parent, "orphan", "orphan", micom_python)
    orphan_medium = orphan / "data" / "glucose-plus-orphan.csv"
    orphan_medium.write_text(
        "nutrient_id,uptake\nglc__D_e,0.5\no2_e,0.5\n", encoding="utf-8"
    )
    orphan_config = _write_config(
        orphan,
        source="predefined",
        aggregation="sum",
        predefined_medium=str(orphan_medium),
        min_growth=0.1,
        open_exchanges=True,
        exchange_bound=1000,
    )
    orphan_source = _source_snapshot(orphan)
    orphan_run = _assert_success(_run(orphan, orphan_config, conda_prefix, "01-orphan"))
    assert "localrule discover_individual_medium:" not in orphan_run
    orphan_results = orphan / "results" / "native"
    medium_status = json.loads(
        (orphan_results / "media" / "custom.status.json").read_text()
    )
    assert medium_status["unused_constituents"] == ["o2_e"]
    assert medium_status["unavailable_models_by_nutrient"] == {
        "o2_e": ["orphan-one", "orphan-two"]
    }
    mapping = _rows(orphan_results / "media" / "custom.mapping.tsv", "\t")
    assert [row["nutrient_id"] for row in mapping] == ["glc__D_e", "glc__D_e"]
    for tool in ("micom", "steadycom"):
        status = json.loads(
            (orphan_results / tool / "custom" / "status.json").read_text()
        )
        assert status["synthetic"] is False and status["threads"] == 1
    micom_growth = _rows(
        orphan_results / "micom" / "custom" / "community_growth.tsv", "\t"
    )
    steady_growth = _rows(
        orphan_results / "steadycom" / "custom" / "community_growth.tsv", "\t"
    )
    assert all(float(row["growth"]) > 0 for row in micom_growth)
    assert float(steady_growth[0]["growth"]) > 0
    micom_ac = [
        float(row["flux"])
        for row in _rows(orphan_results / "micom" / "custom" / "fluxes.tsv", "\t")
        if "EX_ac_e" in row["reaction_id"]
    ]
    steady_ac = [
        float(row["flux"])
        for row in _rows(
            orphan_results / "steadycom" / "custom" / "exchanges.tsv", "\t"
        )
        if "ac_e" in row["reaction_id"] or row["metabolite_id"] == "M_ac_e"
    ]
    assert micom_ac and all(abs(value) <= 1e-9 for value in micom_ac)
    assert steady_ac and all(abs(value) <= 1e-9 for value in steady_ac)
    _assert_source_unchanged(orphan_source)

    unusable = _new_case(native_parent, "unusable", "orphan", micom_python)
    unusable_medium = unusable / "data" / "orphan-only.csv"
    unusable_medium.write_text("nutrient_id,uptake\no2_e,1\n", encoding="utf-8")
    unusable_config = _write_config(
        unusable,
        source="predefined",
        aggregation="sum",
        predefined_medium=str(unusable_medium),
        min_growth=0.1,
        open_exchanges=True,
        exchange_bound=1000,
    )
    unusable_source = _source_snapshot(unusable)
    unusable_result = _run(
        unusable,
        unusable_config,
        conda_prefix,
        "01-unusable",
        "results/native/media/custom.csv",
    )
    assert unusable_result.returncode != 0
    assert "None of the selected medium nutrients maps" in (
        unusable_result.stdout + unusable_result.stderr
    )
    _assert_source_unchanged(unusable_source)

    (main / "native-regression-passed.json").write_text(
        json.dumps(
            {
                "passed": True,
                "conda_prefix": str(conda_prefix),
                "micom_environment": micom_python.parents[1].name,
                "steadycom_environment": steady_python.parents[1].name,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
