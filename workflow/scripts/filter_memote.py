import sys

sys.path.insert(0, str(snakemake.scriptdir))

from common.memote import MISSING_SCORE, memote_model_id, read_memote_report_row

score_files = snakemake.input.scores
quality_out = snakemake.output.quality
filtered_out = snakemake.output.filtered
log_file = snakemake.log[0]

min_total_score = float(snakemake.params.min_total_score)

if isinstance(score_files, str):
    score_files = [score_files]


with open(log_file, "w", encoding="utf-8") as logf:
    logf.write(f"Memote quality filter (min_total_score={min_total_score})\n\n")
    rows = []

    for path in score_files:
        mag = memote_model_id(path)
        try:
            row, section_names = read_memote_report_row(path)
            total = row["total_score"]
            passed = isinstance(total, (int, float)) and total >= min_total_score
            row["status"] = "PASS" if passed else "FAIL"
            rows.append(row)
            logf.write(f"  {mag}: total={total} -> {'PASS' if passed else 'FAIL'}\n")
            logf.write(f"    sections: {list(section_names)}\n")

        except Exception as error:
            logf.write(f"  ERROR reading {path}: {error}\n")
            rows.append(
                {
                    "mag": mag,
                    "total_score": MISSING_SCORE,
                    "consistency": MISSING_SCORE,
                    "annotation": MISSING_SCORE,
                    "status": "ERROR",
                }
            )

    header = ["mag", "total_score", "consistency", "annotation", "status"]
    with open(quality_out, "w", encoding="utf-8") as out:
        out.write("\t".join(header) + "\n")
        for r in rows:
            out.write("\t".join(str(r[h]) for h in header) + "\n")

    passed_models = [r["mag"] for r in rows if r["status"] == "PASS"]
    with open(filtered_out, "w", encoding="utf-8") as out:
        for mag in passed_models:
            out.write(mag + "\n")

    logf.write(f"\nSummary: {len(passed_models)}/{len(rows)} models passed\n")
