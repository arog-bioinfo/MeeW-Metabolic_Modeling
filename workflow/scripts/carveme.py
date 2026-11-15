import importlib.util
import shlex
import sys
from pathlib import Path

from snakemake.shell import shell

sys.path.insert(0, str(snakemake.scriptdir))

from cplex_setup import ensure_cplex_ready
from meew_metabolic_direct_inputs import annotation_faa


def _load_staging_helper():
    source = Path(snakemake.scriptdir) / "common" / "carveme.py"
    spec = importlib.util.spec_from_file_location(
        "_meew_metabolic_carveme_staging", source
    )
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load CarveMe staging helper from {source}")
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    return helper


def _load_index_helper():
    source = Path(snakemake.scriptdir) / "common" / "carveme_index.py"
    spec = importlib.util.spec_from_file_location(
        "_meew_metabolic_carveme_index", source
    )
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load CarveMe index helper from {source}")
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    return helper


carveme_helper = _load_staging_helper()
staged_faa = carveme_helper.staged_faa
carveme_index_helper = _load_index_helper()


def _optional_named_input(name):
    value = getattr(snakemake.input, name, "")
    if not value:
        return ""
    if isinstance(value, str):
        return value
    values = list(value)
    if len(values) != 1:
        raise ValueError(
            f"CarveMe input {name!r} must contain at most one path; found {values!r}."
        )
    return str(values[0])


genome_input = snakemake.input.genome
genome = (
    annotation_faa(genome_input)
    if bool(snakemake.params.get("direct_annotation", False))
    else genome_input
)
model = snakemake.output.model
solver = snakemake.params.get("solver", "scip")
media = snakemake.params.get("media", "M9")
mock = bool(snakemake.params.get("mock", False))
mag = str(snakemake.params.get("mag", getattr(snakemake.wildcards, "mag", ""))).strip()
faa_samples = _optional_named_input("faa_samples")
gtdb_classification = _optional_named_input("gtdb_classification")
index_ready = _optional_named_input("index_ready")

explicit_universe = None
explicit_fallback_reason = ""
if faa_samples:
    sample_records = carveme_helper.read_faa_samples(faa_samples)
    if mag not in sample_records:
        raise ValueError(
            f"FAA sample TSV {faa_samples} has no exact row for requested MAG {mag!r}."
        )
    sample_record = sample_records[mag]
    explicit_universe = sample_record["universe"]
    explicit_fallback_reason = sample_record["universe_fallback_reason"]

universe_selection = carveme_helper.resolve_universe(
    mag,
    explicit_universe=explicit_universe,
    explicit_source=faa_samples,
    explicit_fallback_reason=explicit_fallback_reason,
    gtdb_classification=gtdb_classification,
)
universe = universe_selection.universe

provenance = (
    f"CarveMe universe selection: MAG={mag} universe={universe} "
    f"source={universe_selection.source}"
)
if universe_selection.source_path:
    provenance += f" source_path={universe_selection.source_path}"
if universe_selection.reason:
    provenance += f" fallback_reason={universe_selection.reason}"

if snakemake.log:
    log_path = Path(snakemake.log[0])
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text(provenance + "\n", encoding="utf-8")
    log = f">> {shlex.quote(str(log_path))} 2>&1"
else:
    log = ""

if mock:
    Path(model).parent.mkdir(parents=True, exist_ok=True)
    with open(model, "w") as handle:
        handle.write(
            """<?xml version='1.0' encoding='UTF-8'?><sbml level='3' version='1'></sbml>\n"""
        )
    if snakemake.log:
        with open(snakemake.log[0], "a", encoding="utf-8") as handle:
            handle.write(f"Mock CarveMe output for {genome}\n")
    raise SystemExit

if not index_ready:
    raise RuntimeError(
        "Native CarveMe execution requires this workflow invocation's verified "
        "index readiness, but no initialization output was declared."
    )
carveme_index_helper.verify_readiness(index_ready)

cplex_config = snakemake.config.get("cplex", {})
ensure_cplex_ready(
    solver,
    cplex_config.get("studio_dir", ""),
    auto_upgrade=bool(cplex_config.get("auto_upgrade", True)),
    log_file=snakemake.log[0] if snakemake.log else None,
)

with staged_faa(genome, snakemake.resources.tmpdir) as staged_genome:
    shell(
        "carve {staged_genome:q} "
        "--output {model:q} "
        "--fbc2 "
        "--universe {universe:q} "
        "--solver {solver:q} "
        "--gapfill {media:q} "
        "{log}"
    )

model_path = Path(model)
if not model_path.is_file() or model_path.stat().st_size == 0:
    log_reference = str(snakemake.log[0]) if snakemake.log else "(no log configured)"
    raise RuntimeError(
        "CarveMe returned without producing a non-empty model at "
        f"{model_path}. See the CarveMe log: {log_reference}"
    )
