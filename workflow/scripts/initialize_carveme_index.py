import importlib.util
from pathlib import Path


def _load_helper(path):
    source = Path(path)
    spec = importlib.util.spec_from_file_location(
        "_meew_metabolic_carveme_index", source
    )
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load CarveMe index helper from {source}")
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    return helper


helper = _load_helper(snakemake.input.helper)
log_path = snakemake.log[0] if snakemake.log else None
paths = helper.discover_bundled_index()
status, config_identity, fasta_identity, database_identity = (
    helper.ensure_bundled_index(paths, log_path=log_path)
)
helper.write_readiness(
    snakemake.output[0],
    paths=paths,
    status=status,
    config_identity=config_identity,
    fasta_identity=fasta_identity,
    database_identity=database_identity,
    initializer_log=log_path,
)
