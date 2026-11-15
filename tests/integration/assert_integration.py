#!/usr/bin/env python3
"""Assert contracts for the real-tool integration scenarios."""

import argparse
import csv
import json
import math
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

FIXTURE_IDS = (
    "ecoli_k12_mg1655_GCF_000005845.2",
    "ecoli_k12_w3110_GCF_000010245.2",
)
FILTER_STATUSES = {
    "ecoli_k12_mg1655_GCF_000005845.2": "PASS",
    "ecoli_k12_w3110_GCF_000010245.2": "FAIL",
}
ALL_PASS_STATUSES = {fixture_id: "PASS" for fixture_id in FIXTURE_IDS}
MOCK_PATTERN = re.compile(r"\bmock\b", re.IGNORECASE)


class Validation:
    def __init__(self):
        self.failures = []
        self.passes = []
        self.skips = []
        self.notes = []

    def check(self, label, function, *args):
        try:
            detail = function(*args)
        except (
            AssertionError,
            KeyError,
            OSError,
            TypeError,
            ValueError,
            ET.ParseError,
        ) as error:
            self.failures.append(f"{label}: {error}")
        else:
            self.passes.append(f"{label}: {detail}")

    def require_files(self, label, paths, require_complete):
        existing = [path for path in paths if path.is_file()]
        if len(existing) == len(paths):
            return True
        missing = [str(path) for path in paths if not path.is_file()]
        detail = f"{len(existing)}/{len(paths)} artifacts present; missing {missing}"
        if require_complete:
            self.failures.append(f"{label}: {detail}")
        else:
            self.skips.append(f"{label}: {detail}")
        return False

    def report(self):
        for detail in self.passes:
            print(f"PASS: {detail}")
        for detail in self.notes:
            print(f"NOTE: {detail}")
        for detail in self.skips:
            print(f"SKIP: {detail}")
        for detail in self.failures:
            print(f"FAIL: {detail}", file=sys.stderr)
        print(
            "SUMMARY: "
            f"{len(self.passes)} passed, {len(self.failures)} failed, "
            f"{len(self.skips)} skipped, {len(self.notes)} notes"
        )
        return 1 if self.failures else 0


def assert_regular_file(path):
    assert path.is_file(), f"missing file: {path}"
    assert path.stat().st_size > 0, f"empty file: {path}"


def reject_mock_text(path):
    text = path.read_text(encoding="utf-8", errors="replace")
    assert not MOCK_PATTERN.search(text), f"mock/placeholder marker found in {path}"
    return text


def read_tsv(path, expected_header=None):
    assert_regular_file(path)
    text = reject_mock_text(path)
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = reader.fieldnames
        assert header is not None and len(header) >= 2, f"invalid TSV header in {path}"
        if expected_header is not None:
            assert (
                header == expected_header
            ), f"unexpected header in {path}: expected {expected_header}, found {header}"
        rows = list(reader)
    assert "\t" in text.splitlines()[0], f"header is not tab-delimited in {path}"
    assert rows, f"TSV has no data rows: {path}"
    assert all(None not in row for row in rows), f"malformed TSV row in {path}"
    assert all(
        any(value not in (None, "") for value in row.values()) for row in rows
    ), f"empty TSV row in {path}"
    return header, rows


def assert_exact_ids(rows, id_column, expected_ids=FIXTURE_IDS):
    observed = [row.get(id_column, "") for row in rows]
    assert len(observed) == len(set(observed)), f"duplicate IDs: {observed}"
    assert set(observed) == set(
        expected_ids
    ), f"expected IDs {sorted(expected_ids)}, found {sorted(observed)}"


def assert_qa(results, expected_statuses):
    summary = results / "qa" / "external_qa_summary.tsv"
    filtered = results / "qa" / "filtered_mags.txt"
    _, rows = read_tsv(
        summary,
        ["mag", "faa", "completeness", "contamination", "status"],
    )
    assert_exact_ids(rows, "mag")
    observed_statuses = {row["mag"]: row["status"] for row in rows}
    assert (
        observed_statuses == expected_statuses
    ), f"expected statuses {expected_statuses}, found {observed_statuses}"
    for row in rows:
        assert (
            Path(row["faa"]).stem == row["mag"]
        ), f"FAA stem does not match fixture ID for {row['mag']}: {row['faa']}"
        for column in ("completeness", "contamination"):
            value = float(row[column])
            assert math.isfinite(value), f"non-finite {column} for {row['mag']}"

    assert_regular_file(filtered)
    passed_ids = [
        line.strip() for line in filtered.read_text(encoding="utf-8").splitlines()
    ]
    expected_passed = sorted(
        fixture_id
        for fixture_id, status in expected_statuses.items()
        if status == "PASS"
    )
    assert (
        passed_ids == expected_passed
    ), f"expected filtered IDs {expected_passed}, found {passed_ids}"
    return f"exact statuses and filtered IDs verified: {observed_statuses}"


def local_name(element):
    return element.tag.rsplit("}", 1)[-1]


def assert_model(path):
    assert_regular_file(path)
    reject_mock_text(path)
    root = ET.parse(path).getroot()
    assert local_name(root) == "sbml", f"XML root is not SBML in {path}"
    models = [element for element in root.iter() if local_name(element) == "model"]
    reactions = [
        element for element in root.iter() if local_name(element) == "reaction"
    ]
    assert models, f"SBML has no model element: {path}"
    assert reactions, f"SBML has no reactions and appears to be a placeholder: {path}"
    return f"nonempty parseable SBML with {len(reactions)} reactions"


def numeric(value, label):
    assert not isinstance(value, bool), f"{label} is boolean, not numeric"
    number = float(value)
    assert math.isfinite(number), f"{label} is not finite: {value}"
    return number


def assert_memote_pair(report, score, fixture_id):
    assert_regular_file(report)
    report_text = reject_mock_text(report)
    assert "<html" in report_text.lower(), f"Memote report is not HTML: {report}"
    assert_regular_file(score)
    reject_mock_text(score)
    with score.open(encoding="utf-8") as handle:
        payload = json.load(handle)
    total_score = numeric(
        payload.get("score", {}).get("total_score"),
        f"Memote total score for {fixture_id}",
    )
    return f"real report and numeric total score {total_score}"


def assert_memote_summary(path):
    _, rows = read_tsv(path)
    assert_exact_ids(rows, "mag")
    scores = {row["mag"]: numeric(row["total_score"], row["mag"]) for row in rows}
    return f"numeric total scores for both fixture IDs: {scores}"


def assert_memote_quality(path):
    _, rows = read_tsv(path)
    assert_exact_ids(rows, "mag")
    scores = {row["mag"]: numeric(row["total_score"], row["mag"]) for row in rows}
    statuses = {row["mag"]: row["status"] for row in rows}
    assert (
        statuses == ALL_PASS_STATUSES
    ), f"expected Memote statuses {ALL_PASS_STATUSES}, found {statuses}"
    return f"numeric scores and all-pass statuses verified: {scores}"


def assert_filtered_models(path):
    assert_regular_file(path)
    observed = [line.strip() for line in path.read_text(encoding="utf-8").splitlines()]
    assert observed == sorted(
        FIXTURE_IDS
    ), f"expected both Memote-passing IDs, found {observed}"
    return f"both expected fixture IDs: {observed}"


def assert_generic_tsv(path):
    header, rows = read_tsv(path)
    return f"parseable nonempty TSV ({len(rows)} rows; columns {header})"


def assert_selected_medium(medium_path, mapping_path, status_path, tag):
    assert_regular_file(medium_path)
    with medium_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        assert reader.fieldnames == ["nutrient_id", "uptake"]
        medium_rows = list(reader)
    assert medium_rows and all(None not in row for row in medium_rows)
    nutrients = {row["nutrient_id"] for row in medium_rows}
    assert all(numeric(row["uptake"], row["nutrient_id"]) > 0 for row in medium_rows)

    _, mapping_rows = read_tsv(
        mapping_path,
        [
            "mag",
            "nutrient_id",
            "raw_species_id",
            "cobra_exchange_id",
            "exchange_coefficient",
        ],
    )
    assert {row["mag"] for row in mapping_rows} == set(FIXTURE_IDS)
    assert {row["nutrient_id"] for row in mapping_rows} <= nutrients
    assert all(
        math.isclose(abs(numeric(row["exchange_coefficient"], "coefficient")), 1.0)
        for row in mapping_rows
    )

    assert_regular_file(status_path)
    status = json.loads(status_path.read_text(encoding="utf-8"))
    assert status["synthetic"] is False
    assert status["tag"] == tag and status["source"] == "discover"
    assert status["model_ids"] == list(FIXTURE_IDS)
    return f"{len(nutrients)}-nutrient native {tag} medium and raw-source mapping"


def assert_micom_growth(path, require_complete):
    _, rows = read_tsv(
        path,
        ["stage", "growth", "objective_value", "stage_status", "synthetic"],
    )
    assert {row["stage"] for row in rows} == {"fba", "pfba"}
    assert all(row["synthetic"] == "False" for row in rows)
    partial = False
    for row in rows:
        if row["stage_status"] == "complete":
            numeric(row["growth"], f"MICOM {row['stage']} growth")
            numeric(row["objective_value"], f"MICOM {row['stage']} objective")
        else:
            assert row["stage"] == "pfba"
            assert row["stage_status"] == "unavailable"
            assert row["growth"] == row["objective_value"] == "NA"
            partial = True
    if partial:
        assert not require_complete, "MICOM pFBA is explicitly unavailable"
        return "validated native MICOM FBA with unavailable pFBA marked NA"
    return "finite native MICOM FBA/pFBA community growth"


def assert_micom_members(path, require_complete):
    _, rows = read_tsv(path, ["stage", "mag", "abundance", "growth", "synthetic"])
    stages = {row["stage"] for row in rows}
    assert stages in ({"fba"}, {"fba", "pfba"})
    expected = {(stage, fixture_id) for stage in stages for fixture_id in FIXTURE_IDS}
    assert {(row["stage"], row["mag"]) for row in rows} == expected
    for row in rows:
        assert row["synthetic"] == "False"
        assert numeric(row["abundance"], f"MICOM abundance for {row['mag']}") >= 0
        numeric(row["growth"], f"MICOM growth for {row['mag']}")
    if stages == {"fba"}:
        assert not require_complete, "MICOM pFBA member results are unavailable"
        return "complete native MICOM FBA members; pFBA unavailable"
    return "complete native MICOM FBA/pFBA member growth"


def assert_micom_fluxes(path, require_complete):
    _, rows = read_tsv(
        path, ["stage", "compartment", "reaction_id", "flux", "synthetic"]
    )
    stages = {row["stage"] for row in rows}
    assert stages in ({"fba"}, {"fba", "pfba"})
    assert all(row["synthetic"] == "False" for row in rows)
    fluxes = [numeric(row["flux"], f"MICOM flux {row['reaction_id']}") for row in rows]
    assert any(value != 0 for value in fluxes), "MICOM returned only zero fluxes"
    if stages == {"fba"}:
        assert not require_complete, "MICOM pFBA fluxes are unavailable"
        return f"{len(rows)} finite native FBA flux values; pFBA unavailable"
    return f"{len(rows)} finite native MICOM flux values"


def assert_steadycom_growth(path):
    _, rows = read_tsv(path, ["stage", "growth", "synthetic"])
    assert len(rows) == 1 and rows[0]["stage"] == "steadycom"
    assert rows[0]["synthetic"] == "False"
    assert numeric(rows[0]["growth"], "SteadyCom growth") > 0
    return "finite positive native SteadyCom growth"


def assert_steadycom_abundances(path):
    _, rows = read_tsv(path, ["mag", "abundance", "synthetic"])
    assert_exact_ids(rows, "mag")
    assert all(row["synthetic"] == "False" for row in rows)
    abundances = {
        row["mag"]: numeric(row["abundance"], f"SteadyCom abundance for {row['mag']}")
        for row in rows
    }
    assert all(
        value >= 0 for value in abundances.values()
    ), f"negative SteadyCom abundance: {abundances}"
    assert math.isclose(sum(abundances.values()), 1.0, abs_tol=1e-6)
    return f"numeric abundances for expected fixture IDs: {abundances}"


def assert_steadycom_exchanges(path):
    _, rows = read_tsv(
        path,
        ["scope", "mag", "metabolite_id", "reaction_id", "flux", "synthetic"],
    )
    assert {row["scope"] for row in rows} <= {"environment", "member"}
    assert all(row["synthetic"] == "False" for row in rows)
    fluxes = [numeric(row["flux"], "SteadyCom exchange flux") for row in rows]
    assert any(value != 0 for value in fluxes), "SteadyCom returned only zero exchanges"
    return f"{len(rows)} finite native SteadyCom exchange values"


def assert_steadycom_variability(path, require_complete):
    _, rows = read_tsv(
        path,
        [
            "mag",
            "minimum_abundance",
            "maximum_abundance",
            "objective_fraction",
            "variability_status",
            "synthetic",
        ],
    )
    assert_exact_ids(rows, "mag")
    missing = []
    unavailable = False
    for row in rows:
        bounds = {}
        for direction, column in (
            ("minimum", "minimum_abundance"),
            ("maximum", "maximum_abundance"),
        ):
            if row[column] == "NA":
                bounds[direction] = None
                missing.append((row["mag"], direction))
            else:
                bounds[direction] = numeric(
                    row[column], f"SteadyComVA {direction} bound"
                )
                assert 0 <= bounds[direction] <= 1
        if all(value is not None for value in bounds.values()):
            assert bounds["minimum"] <= bounds["maximum"]
            assert row["variability_status"] == "complete"
        elif all(value is None for value in bounds.values()):
            assert row["variability_status"] in {"partial", "unavailable"}
            unavailable = unavailable or row["variability_status"] == "unavailable"
        else:
            assert row["variability_status"] == "partial"
        assert math.isclose(
            numeric(row["objective_fraction"], "objective fraction"), 0.9
        )
        assert row["synthetic"] == "False"
    if unavailable:
        assert all(row["variability_status"] == "unavailable" for row in rows)
        assert not require_complete, "SteadyComVA is explicitly unavailable"
        return "explicit unavailable native SteadyComVA output with NA bounds"
    if missing:
        assert not require_complete, f"SteadyComVA is partial; missing {missing}"
        return f"explicit partial native SteadyComVA output; missing {missing}"
    return "finite native SteadyComVA 0.9 intervals for both fixture IDs"


def assert_native_status(path, tool, tag):
    assert_regular_file(path)
    payload = json.loads(path.read_text(encoding="utf-8"))
    assert payload["synthetic"] is False
    assert payload["tool"] == tool and payload["medium_tag"] == tag
    assert payload["model_ids"] == list(FIXTURE_IDS)
    return f"native {tool} status for {tag} and the complete cohort"


def assert_micom_status(path, tag, require_complete):
    assert_regular_file(path)
    payload = json.loads(path.read_text(encoding="utf-8"))
    assert payload["synthetic"] is False
    assert payload["tool"] == "micom" and payload["medium_tag"] == tag
    assert payload["model_ids"] == list(FIXTURE_IDS)
    assert payload["analysis_status"] in {"complete", "partial"}
    stages = payload["stage_status"]
    assert stages["fba"]["result_status"] == "complete"
    assert stages["pfba"]["result_status"] in {"complete", "unavailable"}
    for stage in ("fba", "pfba"):
        assert "native_status" in stages[stage]
        assert "native_message" in stages[stage]
        assert "exception_type" in stages[stage]
        assert "exception_message" in stages[stage]
    if stages["pfba"]["result_status"] == "unavailable":
        assert payload["analysis_status"] == "partial"
        assert not require_complete, "MICOM pFBA is explicitly unavailable"
        return f"native MICOM status for {tag} with unavailable pFBA"
    assert payload["analysis_status"] == "complete"
    return f"native MICOM status for {tag} with complete FBA/pFBA"


def assert_steadycom_status(path, tag, require_complete):
    assert_regular_file(path)
    payload = json.loads(path.read_text(encoding="utf-8"))
    assert payload["synthetic"] is False
    assert payload["tool"] == "steadycom" and payload["medium_tag"] == tag
    assert payload["model_ids"] == list(FIXTURE_IDS)
    variability_status = payload["variability_status"]
    missing = payload["missing_variability_bounds"]
    assert payload["analysis_status"] in {"complete", "partial"}
    assert variability_status in {"complete", "partial", "unavailable"}
    if variability_status == "complete":
        assert missing == []
        assert payload["analysis_status"] == "complete"
        return f"native steadycom status for {tag} with complete variability"
    if variability_status == "unavailable":
        assert missing == []
        initialization = payload["stage_status"]["steadycom_va"]["initialization"]
        assert initialization["result_status"] == "unavailable"
        assert "native_status" in initialization
        assert "native_message" in initialization
        assert "exception_type" in initialization
        assert "exception_message" in initialization
        assert payload["analysis_status"] == "partial"
        assert not require_complete, "SteadyComVA is explicitly unavailable"
        return f"native steadycom status for {tag} with unavailable variability"
    assert missing, "partial SteadyComVA status has no missing-bound diagnostics"
    for item in missing:
        assert item["mag"] in FIXTURE_IDS
        assert item["direction"] in {"minimum", "maximum"}
        assert item["bound_status"] in {
            "infeasible",
            "nonoptimal",
            "unknown",
            "unavailable",
        }
        assert "native_status" in item and "native_message" in item
    assert not require_complete, "SteadyComVA is explicitly partial"
    return f"native steadycom status for {tag} with explicit partial variability"


def validate_filter(results, validation):
    validation.check("external QA filter", assert_qa, results, FILTER_STATUSES)


def validate_full(results, validation, require_complete):
    validation.check("external QA all-pass", assert_qa, results, ALL_PASS_STATUSES)

    models = [results / "gems" / f"{fixture_id}.xml" for fixture_id in FIXTURE_IDS]
    validation.require_files("CarveMe model stage", models, require_complete)
    for path in models:
        if path.is_file():
            validation.check(f"CarveMe model {path.stem}", assert_model, path)

    memote_pairs = [
        (
            results / "memote" / fixture_id / "report.html",
            results / "memote" / fixture_id / "score.json",
            fixture_id,
        )
        for fixture_id in FIXTURE_IDS
    ]
    memote_files = [path for pair in memote_pairs for path in pair[:2]]
    validation.require_files("Memote per-model stage", memote_files, require_complete)
    for report, score, fixture_id in memote_pairs:
        if report.is_file() and score.is_file():
            validation.check(
                f"Memote {fixture_id}",
                assert_memote_pair,
                report,
                score,
                fixture_id,
            )

    memote_summary = results / "memote" / "memote_summary.tsv"
    memote_quality = results / "memote" / "memote_quality.tsv"
    filtered_models = results / "memote" / "filtered_models.txt"
    aggregate_files = [memote_summary, memote_quality, filtered_models]
    validation.require_files(
        "Memote aggregate/filter stage", aggregate_files, require_complete
    )
    if memote_summary.is_file():
        validation.check("Memote summary", assert_memote_summary, memote_summary)
    if memote_quality.is_file():
        validation.check("Memote quality filter", assert_memote_quality, memote_quality)
    if filtered_models.is_file():
        validation.check("Memote passing IDs", assert_filtered_models, filtered_models)

    smetana_global = results / "smetana" / "global.tsv"
    smetana_detailed = results / "smetana" / "detailed.tsv"
    smetana_files = [smetana_global, smetana_detailed]
    validation.require_files("SMETANA stage", smetana_files, require_complete)
    for path in smetana_files:
        if path.is_file():
            validation.check(f"SMETANA {path.name}", assert_generic_tsv, path)

    medium_tag = "summed_demand"
    medium_root = results / "media"
    medium_files = [
        medium_root / f"{medium_tag}.csv",
        medium_root / f"{medium_tag}.mapping.tsv",
        medium_root / f"{medium_tag}.status.json",
    ]
    validation.require_files("selected medium stage", medium_files, require_complete)
    if all(path.is_file() for path in medium_files):
        validation.check(
            "selected simulation medium",
            assert_selected_medium,
            *medium_files,
            medium_tag,
        )

    micom_root = results / "micom" / medium_tag
    micom_files = {
        "community growth": micom_root / "community_growth.tsv",
        "member growth": micom_root / "member_growth.tsv",
        "fluxes": micom_root / "fluxes.tsv",
        "status": micom_root / "status.json",
    }
    validation.require_files(
        "MICOM stage", list(micom_files.values()), require_complete
    )
    if micom_files["community growth"].is_file():
        validation.check(
            "MICOM community growth",
            assert_micom_growth,
            micom_files["community growth"],
            require_complete,
        )
    if micom_files["member growth"].is_file():
        validation.check(
            "MICOM member growth",
            assert_micom_members,
            micom_files["member growth"],
            require_complete,
        )
    if micom_files["fluxes"].is_file():
        validation.check(
            "MICOM fluxes", assert_micom_fluxes, micom_files["fluxes"], require_complete
        )
    if micom_files["status"].is_file():
        validation.check(
            "MICOM status",
            assert_micom_status,
            micom_files["status"],
            medium_tag,
            require_complete,
        )

    steadycom_root = results / "steadycom" / medium_tag
    steadycom_files = {
        "community growth": steadycom_root / "community_growth.tsv",
        "abundances": steadycom_root / "abundances.tsv",
        "exchanges": steadycom_root / "exchanges.tsv",
        "variability": steadycom_root / "variability.tsv",
        "status": steadycom_root / "status.json",
    }
    validation.require_files(
        "SteadyCom stage", list(steadycom_files.values()), require_complete
    )
    if steadycom_files["community growth"].is_file():
        validation.check(
            "SteadyCom community growth",
            assert_steadycom_growth,
            steadycom_files["community growth"],
        )
    if steadycom_files["abundances"].is_file():
        validation.check(
            "SteadyCom abundances",
            assert_steadycom_abundances,
            steadycom_files["abundances"],
        )
    if steadycom_files["exchanges"].is_file():
        validation.check(
            "SteadyCom exchanges",
            assert_steadycom_exchanges,
            steadycom_files["exchanges"],
        )
    if steadycom_files["variability"].is_file():
        validation.check(
            "SteadyCom variability",
            assert_steadycom_variability,
            steadycom_files["variability"],
            require_complete,
        )
    if steadycom_files["status"].is_file():
        validation.check(
            "SteadyCom status",
            assert_steadycom_status,
            steadycom_files["status"],
            medium_tag,
            require_complete,
        )


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("scenario", choices=("filter", "full"))
    parser.add_argument("--results", required=True, type=Path)
    parser.add_argument(
        "--require-complete",
        action="store_true",
        help="Fail if any full-scenario stage is absent.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    validation = Validation()
    if args.scenario == "filter":
        if args.require_complete:
            validation.failures.append(
                "--require-complete is only meaningful for the full scenario"
            )
        validate_filter(args.results.resolve(), validation)
    else:
        validate_full(args.results.resolve(), validation, args.require_complete)
    return validation.report()


if __name__ == "__main__":
    raise SystemExit(main())
