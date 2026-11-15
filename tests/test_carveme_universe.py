import os
import subprocess
import sys
import time
from pathlib import Path

import pytest
import yaml

SCRIPTS = Path(__file__).resolve().parents[1] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

from common.carveme import (
    read_faa_samples,
    read_gtdb_domains,
    resolve_universe,
)

FIXTURE = Path(__file__).parent / "fixtures" / "gtdbtk_classification.tsv"
ROOT = Path(__file__).resolve().parents[1]
SNAKEFILE = ROOT / "workflow" / "Snakefile"


def test_two_column_faa_sheet_defaults_and_preserves_expanded_relative_paths(
    tmp_path, monkeypatch
):
    monkeypatch.setenv("FAA_ROOT", "relative inputs")
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "mag\tpath\n" "  mag-b  \t$FAA_ROOT/b.faa\n" "mag-a\t~/proteins/a.faa\n",
        encoding="utf-8",
    )

    records = read_faa_samples(sheet)

    assert list(records) == ["mag-b", "mag-a"]
    assert records["mag-b"] == {
        "path": "relative inputs/b.faa",
        "universe": None,
        "universe_fallback_reason": "FAA sample TSV has no universe column",
    }
    assert records["mag-a"]["path"] == str(Path.home() / "proteins/a.faa")
    selection = resolve_universe(
        "mag-b",
        explicit_source=str(sheet),
        explicit_fallback_reason=records["mag-b"]["universe_fallback_reason"],
    )
    assert selection.universe == "bacteria"
    assert selection.source == "default"
    assert selection.source_path == str(sheet)


def test_three_column_faa_sheet_supports_mixed_and_blank_universes(tmp_path):
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "mag\tpath\tuniverse\tnote\n"
        "mag-b\tb.faa\t bacteria \tbacterial\n"
        "mag-a\ta.faa\tarchaea\tarchaeal\n"
        "mag-default\td.faa\t   \tblank\n",
        encoding="utf-8",
    )

    records = read_faa_samples(sheet)

    assert {mag: record["universe"] for mag, record in records.items()} == {
        "mag-b": "bacteria",
        "mag-a": "archaea",
        "mag-default": None,
    }
    assert (
        records["mag-default"]["universe_fallback_reason"]
        == "blank universe cell in FAA sample TSV"
    )


@pytest.mark.parametrize("value", ["archaeal", "Bacteria", "virus"])
def test_invalid_nonblank_sheet_universe_reports_mag_and_line(tmp_path, value):
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        f"mag\tpath\tuniverse\nmag-invalid\tinput.faa\t{value}\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match=r"mag-invalid.*line 2.*bacteria.*archaea"):
        read_faa_samples(sheet)


@pytest.mark.parametrize(
    ("contents", "message"),
    [
        ("mag\tuniverse\nmag-a\tbacteria\n", "missing required column.*path"),
        ("mag\tpath\nmag-a\ta.faa\nmag-a\tb.faa\n", "duplicate mag.*line 3"),
        ("mag\tpath\nmag-a\n", "malformed row"),
        ("mag\tpath\n\ta.faa\n", "empty mag/path.*line 2"),
    ],
)
def test_faa_sheet_preserves_strict_identity_and_structure_validation(
    tmp_path, contents, message
):
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(contents, encoding="utf-8")
    with pytest.raises(ValueError, match=message):
        read_faa_samples(sheet)


@pytest.mark.parametrize(
    ("mag", "expected"),
    [
        ("EST6_prok_bin_002", "bacteria"),
        ("EST6_prok_bin_005", "archaea"),
        ("EST6_prok_bin_010", "archaea"),
    ],
)
def test_est6_shaped_gtdb_fixture_maps_exact_domains(mag, expected):
    selection = resolve_universe(mag, gtdb_classification=FIXTURE)
    assert selection.universe == expected
    assert selection.source == "GTDB"
    assert selection.source_path == str(FIXTURE)
    assert selection.reason == ""


@pytest.mark.parametrize(
    ("mag", "reason"),
    [
        ("EST6_prok_bin_unknown", "blank GTDB-Tk domain"),
        ("EST6_prok_bin_other", "unsupported GTDB-Tk domain 'unknown'"),
        ("EST6_prok_bin_005.faa", "no exact GTDB-Tk user_genome row"),
        ("not-present", "no exact GTDB-Tk user_genome row"),
    ],
)
def test_gtdb_unknown_blank_and_nonmatching_ids_fall_back_visibly(mag, reason):
    selection = resolve_universe(mag, gtdb_classification=FIXTURE)
    assert selection.universe == "bacteria"
    assert selection.source == "GTDB fallback"
    assert reason in selection.reason


def test_explicit_choice_has_priority_without_reading_taxonomy(tmp_path):
    missing = tmp_path / "future-classification.tsv"
    selection = resolve_universe(
        "mag-a",
        explicit_universe="archaea",
        explicit_source="samples.tsv",
        gtdb_classification=missing,
    )
    assert selection == ("archaea", "explicit TSV", "samples.tsv", "")


@pytest.mark.parametrize(
    ("contents", "message"),
    [
        ("user_genome\nmag-a\n", "missing required column.*domain"),
        ("domain\tuser_genome\nbacteria\n", "malformed row"),
        ("domain\tuser_genome\nbacteria\t\n", "blank user_genome"),
        (
            "domain\tuser_genome\nbacteria\tmag-a\nbacteria\tmag-a\n",
            "duplicate rows.*mag-a",
        ),
        (
            "domain\tuser_genome\nbacteria\tmag-a\narchaea\tmag-a\n",
            "conflicting rows.*mag-a",
        ),
    ],
)
def test_gtdb_table_rejects_invalid_headers_rows_and_duplicate_ids(
    tmp_path, contents, message
):
    table = tmp_path / "classification.tsv"
    table.write_text(contents, encoding="utf-8")
    with pytest.raises(ValueError, match=message):
        read_gtdb_domains(table)


def test_requested_missing_gtdb_table_does_not_fall_back(tmp_path):
    missing = tmp_path / "missing.tsv"
    with pytest.raises(FileNotFoundError, match="(?i)requested.*missing or unreadable"):
        resolve_universe("mag-a", gtdb_classification=missing)


def test_offline_standalone_two_then_three_column_sheet_reruns_mixed_mags(
    tmp_path,
):
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    sources = []
    for mag in ("mag-b", "mag-a"):
        source = inputs / f"{mag}.faa"
        source.write_text(f">{mag}\nMPEPTIDE\n", encoding="utf-8")
        sources.append(source)
    before = {
        source: (source.read_bytes(), source.stat().st_mtime_ns) for source in sources
    }
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "mag\tpath\nmag-b\tinputs/mag-b.faa\nmag-a\tinputs/mag-a.faa\n",
        encoding="utf-8",
    )
    config = yaml.safe_load((ROOT / "config" / "config.yaml").read_text())
    config["sample"] = "sheet-sample"
    config["input"] = {
        "faa_samples": str(sheet),
        "faa_globs": [],
        "faa_dir": "",
        "bakta_dir": "",
        "bakta_done": "",
    }
    config["solver"] = "scip"
    config["execution"]["mock"] = True
    config["qa"]["enabled"] = False
    config["community"] = {
        "enabled": False,
        "tools": {"smetana": False, "micom": False, "steadycom": False},
    }
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    targets = [f"results/sheet-sample/gems/{mag}.xml" for mag in ("mag-b", "mag-a")]
    environment = os.environ.copy()
    environment.update(
        CONDA_OFFLINE="true",
        PIP_NO_INDEX="1",
        PYTHONDONTWRITEBYTECODE="1",
        PYTHONPYCACHEPREFIX=str(tmp_path / "pycache"),
        XDG_CACHE_HOME=str(tmp_path / "xdg-cache"),
    )

    def run_workflow():
        return subprocess.run(
            [
                sys.executable,
                "-m",
                "snakemake",
                "--snakefile",
                str(SNAKEFILE),
                "--directory",
                str(tmp_path),
                "--configfile",
                str(config_path),
                "--nolock",
                "--cores",
                "1",
                *targets,
            ],
            cwd=ROOT,
            check=False,
            capture_output=True,
            text=True,
            env=environment,
        )

    first = run_workflow()
    assert first.returncode == 0, first.stdout + first.stderr
    logs = tmp_path / "results" / "sheet-sample" / "logs" / "carveme"
    for mag in ("mag-b", "mag-a"):
        assert f"MAG={mag} universe=bacteria source=default" in (
            logs / f"{mag}.log"
        ).read_text(encoding="utf-8")
    first_model_mtimes = {
        target: (tmp_path / target).stat().st_mtime_ns for target in targets
    }

    time.sleep(0.01)
    sheet.write_text(
        "mag\tpath\tuniverse\n"
        "mag-b\tinputs/mag-b.faa\tbacteria\n"
        "mag-a\tinputs/mag-a.faa\tarchaea\n",
        encoding="utf-8",
    )
    second = run_workflow()
    assert second.returncode == 0, second.stdout + second.stderr
    assert "MAG=mag-b universe=bacteria source=explicit TSV" in (
        logs / "mag-b.log"
    ).read_text(encoding="utf-8")
    assert "MAG=mag-a universe=archaea source=explicit TSV" in (
        logs / "mag-a.log"
    ).read_text(encoding="utf-8")
    assert all(
        (tmp_path / target).stat().st_mtime_ns > first_model_mtimes[target]
        for target in targets
    )
    assert {
        source: (source.read_bytes(), source.stat().st_mtime_ns) for source in sources
    } == before
