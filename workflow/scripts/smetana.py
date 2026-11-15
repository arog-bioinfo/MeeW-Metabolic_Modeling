import sys
from pathlib import Path

sys.path.insert(0, str(snakemake.scriptdir))

from common.community import select_memote_passing_models
from common.smetana import run_smetana_mode, validate_media_database
from cplex_setup import ensure_cplex_ready

models = snakemake.input.models
output = str(snakemake.output[0])
mode = str(snakemake.params.mode)
solver = snakemake.params.get("solver", "scip")
mock = bool(snakemake.params.get("mock", False))
media = list(snakemake.params.get("media", []))
media_db = str(snakemake.input.get("media_db", "") or "")
log = str(snakemake.log[0])

model_paths = select_memote_passing_models(models, snakemake.input.filtered_models)

if not model_paths:
    raise ValueError("No Memote-passing models are available for SMETANA.")

if media_db:
    validate_media_database(media_db, media)

Path(output).parent.mkdir(parents=True, exist_ok=True)
Path(log).parent.mkdir(parents=True, exist_ok=True)

if mock:
    if mode == "global":
        contents = "community\tscore\nmock\t1.0\n"
    elif mode == "detailed":
        contents = "donor\treceiver\tscore\nmock\tmock\t1.0\n"
    else:
        raise ValueError(
            f"Unsupported SMETANA mode {mode!r}; expected 'global' or 'detailed'."
        )
    Path(output).write_text(contents, encoding="utf-8")
    Path(log).write_text(
        f"Mock SMETANA {mode} output from Memote-passing models.\n",
        encoding="utf-8",
    )
    raise SystemExit

# Start this job's log before CPLEX preflight; SMETANA output is appended below.
Path(log).write_text("", encoding="utf-8")
cplex_config = snakemake.config.get("cplex", {})
ensure_cplex_ready(
    solver,
    cplex_config.get("studio_dir", ""),
    auto_upgrade=bool(cplex_config.get("auto_upgrade", True)),
    log_file=log,
)

smetana_path = str(Path(sys.executable).parent / "smetana")
run_smetana_mode(
    executable=smetana_path,
    models=model_paths,
    mode=mode,
    solver=solver,
    output=output,
    log=log,
    tmpdir=snakemake.resources.tmpdir,
    media_db=media_db,
    media=media,
)
