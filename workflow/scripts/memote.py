import json
import os

model_file = snakemake.input.model
report_file = snakemake.output.report
score_file = snakemake.output.score
mock = bool(snakemake.params.get("mock", False))

os.makedirs(os.path.dirname(report_file), exist_ok=True)

if mock:
    with open(report_file, "w", encoding="utf-8") as handle:
        handle.write("<html><body>Mock Memote report</body></html>\n")
    with open(score_file, "w", encoding="utf-8") as handle:
        json.dump(
            {
                "score": {
                    "total_score": 1.0,
                    "sections": [{"section": "consistency", "score": 1.0}],
                }
            },
            handle,
        )
    raise SystemExit

import cobra
from memote.suite.api import test_model
from memote.suite.reporting import ReportConfiguration, SnapshotReport

model = cobra.io.read_sbml_model(model_file)

_, result = test_model(model, results=True)

config = ReportConfiguration.load()
report = SnapshotReport(result=result, configuration=config)

with open(report_file, "w", encoding="utf-8") as f:
    f.write(report.render_html())

with open(score_file, "w", encoding="utf-8") as f:
    f.write(report.render_json())
