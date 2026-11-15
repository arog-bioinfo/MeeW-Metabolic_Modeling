import sys

sys.path.insert(0, str(snakemake.scriptdir))

from common.memote import MISSING_SCORE, memote_model_id, read_memote_report_row

score_files = snakemake.input.scores
summary_out = snakemake.output.summary
log_file = snakemake.log[0]

if isinstance(score_files, str):
    score_files = [score_files]


with open(log_file, "w", encoding="utf-8") as logf:
    rows = []
    for path in score_files:
        mag = memote_model_id(path)
        try:
            row, _ = read_memote_report_row(path)
            row["score_file"] = path
            rows.append(row)
        except Exception as error:
            logf.write(f"ERROR processing {path}: {error}\n")
            rows.append(
                {
                    "mag": mag,
                    "total_score": "ERROR",
                    "consistency": MISSING_SCORE,
                    "annotation": MISSING_SCORE,
                    "score_file": path,
                }
            )

    header = ["mag", "total_score", "consistency", "annotation", "score_file"]
    with open(summary_out, "w", encoding="utf-8") as out:
        out.write("\t".join(header) + "\n")
        for r in rows:
            out.write("\t".join(str(r[h]) for h in header) + "\n")

    logf.write(f"Summary written ({len(rows)} models)\n")
