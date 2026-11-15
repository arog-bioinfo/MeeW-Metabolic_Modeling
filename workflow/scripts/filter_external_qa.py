import csv
import glob
import importlib.util
import math
import os
import sys
from pathlib import Path

sys.path.insert(0, str(snakemake.scriptdir))

from common.path_collections import duplicate_path_stems
from meew_metabolic_direct_inputs import annotation_faa_by_mag


def _load_carveme_input_helper():
    source = Path(snakemake.scriptdir) / "common" / "carveme.py"
    spec = importlib.util.spec_from_file_location(
        "_meew_metabolic_filter_carveme_inputs", source
    )
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load CarveMe input helper from {source}")
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    return helper


read_faa_samples = _load_carveme_input_helper().read_faa_samples

summary_out = snakemake.output.summary
filtered_out = snakemake.output.filtered
log_file = snakemake.log[0]

faa_by_mag = dict(snakemake.params.faa_by_mag)
mags = list(snakemake.params.mags)
bakta_dir = str(snakemake.params.get("bakta_dir", ""))
direct_annotation = bool(snakemake.params.get("direct_annotation", False))
annotation_dirs = snakemake.input.get("annotation_dirs", [])
if isinstance(annotation_dirs, str):
    annotation_dirs = [annotation_dirs]
else:
    annotation_dirs = list(annotation_dirs)
enabled = bool(snakemake.params.enabled)
qa_tsv = str(snakemake.input.get("qa_tsv", "") or "")
faa_samples = str(snakemake.input.get("faa_samples", ""))
id_col = str(snakemake.params.id_col)
completeness_col = str(snakemake.params.completeness_col)
contamination_col = str(snakemake.params.contamination_col)
min_completeness = float(snakemake.params.min_completeness)
max_contamination = float(snakemake.params.max_contamination)

os.makedirs(os.path.dirname(summary_out), exist_ok=True)
os.makedirs(os.path.dirname(log_file), exist_ok=True)


def write_outputs(rows):
    header = ["mag", "faa", "completeness", "contamination", "status"]
    with open(summary_out, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    with open(filtered_out, "w", encoding="utf-8") as handle:
        for row in rows:
            if row["status"] == "PASS":
                handle.write(f"{row['mag']}\n")


def load_qa_rows(path):
    if not path:
        raise ValueError("QA filtering is enabled, but no QA TSV path was provided.")
    required_columns = [id_col, completeness_col, contamination_col]
    qa_rows = {}
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = [
            column
            for column in required_columns
            if column not in (reader.fieldnames or [])
        ]
        if missing:
            raise ValueError(
                f"QA TSV is missing required column(s): {', '.join(missing)}."
            )
        for line_number, record in enumerate(reader, start=2):
            qa_id = record.get(id_col, "")
            if not qa_id:
                raise ValueError(
                    f"QA TSV has an empty ID value in column '{id_col}' on line {line_number}."
                )
            if qa_id in qa_rows:
                raise ValueError(
                    f"QA TSV contains duplicate ID '{qa_id}' on line {line_number}."
                )
            qa_rows[qa_id] = record
    return qa_rows


def load_faa_samples(path):
    return {mag: record["path"] for mag, record in read_faa_samples(path).items()}


def discover_bakta_faa(path):
    files = sorted(
        {
            os.path.normpath(candidate)
            for candidate in glob.glob(
                os.path.join(path, "**", "*.faa"), recursive=True
            )
            if "hypothetical" not in Path(candidate).name.lower()
        }
    )
    duplicates = duplicate_path_stems(files)
    if duplicates:
        raise ValueError(f"Duplicate .faa basenames are not supported: {duplicates}")
    return {Path(faa).stem: faa for faa in files}


if direct_annotation:
    if faa_by_mag or faa_samples or bakta_dir:
        raise ValueError(
            "Direct Annotation directory mode cannot be combined with standalone FAA "
            "sample-sheet or discovery inputs."
        )
    faa_by_mag = annotation_faa_by_mag(mags, annotation_dirs)

if not faa_by_mag and faa_samples:
    faa_by_mag = load_faa_samples(faa_samples)
    mags = sorted(faa_by_mag)

if not faa_by_mag and bakta_dir:
    faa_by_mag = discover_bakta_faa(bakta_dir)
    mags = sorted(faa_by_mag)

if not faa_by_mag:
    raise ValueError(
        "No usable .faa files found. Configure input.faa_samples, input.faa_globs, "
        "input.faa_dir, or input.bakta_dir; Bakta discovery excludes "
        "*hypothetical*.faa files case-insensitively."
    )


with open(log_file, "w", encoding="utf-8") as log:
    if not enabled:
        rows = [
            {
                "mag": mag,
                "faa": faa_by_mag[mag],
                "completeness": "NA",
                "contamination": "NA",
                "status": "PASS",
            }
            for mag in mags
        ]
        write_outputs(rows)
        log.write("External QA filtering is disabled; all discovered FAA files pass.\n")
    else:
        log.write(
            "External QA filtering enabled\n"
            f"  QA TSV: {qa_tsv}\n"
            f"  ID column: {id_col}\n"
            f"  completeness column: {completeness_col}\n"
            f"  contamination column: {contamination_col}\n"
            f"  min completeness: {min_completeness}\n"
            f"  max contamination: {max_contamination}\n"
        )
        qa_rows = load_qa_rows(qa_tsv)
        rows = []
        for mag in mags:
            record = qa_rows.get(mag)
            if record is None:
                rows.append(
                    {
                        "mag": mag,
                        "faa": faa_by_mag[mag],
                        "completeness": "NA",
                        "contamination": "NA",
                        "status": "MISSING",
                    }
                )
                continue
            try:
                completeness = float(record[completeness_col])
                contamination = float(record[contamination_col])
            except (TypeError, ValueError) as error:
                raise ValueError(
                    "External QA TSV contains malformed numeric value(s) "
                    f"for ID '{mag}' in columns '{completeness_col}'/'{contamination_col}'."
                ) from error
            if not (math.isfinite(completeness) and math.isfinite(contamination)):
                raise ValueError(
                    "External QA TSV contains non-finite completeness/contamination "
                    f"for ID '{mag}'."
                )
            passed = (
                completeness >= min_completeness and contamination <= max_contamination
            )
            rows.append(
                {
                    "mag": mag,
                    "faa": faa_by_mag[mag],
                    "completeness": completeness,
                    "contamination": contamination,
                    "status": "PASS" if passed else "FAIL",
                }
            )
        write_outputs(rows)
        passed_count = sum(row["status"] == "PASS" for row in rows)
        log.write(
            f"Summary: {passed_count}/{len(rows)} FAA files passed external QA.\n"
        )
