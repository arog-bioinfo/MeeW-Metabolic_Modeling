"""Symlink discovered FAA inputs into the CheckM2 input directory."""

import importlib.util
import shutil
from pathlib import Path


def _load_direct_inputs_helper():
    source = str(snakemake.input.get("direct_inputs_helper", ""))
    if not source:
        raise ValueError(
            "Bound CheckM2 staging requires its declared direct-input helper."
        )
    spec = importlib.util.spec_from_file_location(
        "_meew_metabolic_stage_direct_inputs", source
    )
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load Metabolic Modeling helper from {source}")
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    return helper


faa_by_mag = dict(snakemake.params.faa_by_mag)
if bool(snakemake.params.get("direct_annotation", False)):
    annotation_faa = _load_direct_inputs_helper().annotation_faa
    faa_by_mag = {
        mag: annotation_faa(Path(faa).parent) for mag, faa in faa_by_mag.items()
    }
if not faa_by_mag and snakemake.params.bakta_dir:
    discovered = sorted(
        path
        for path in Path(snakemake.params.bakta_dir).rglob("*.faa")
        if "hypothetical" not in path.name.lower()
    )
    faa_by_mag = {path.stem: str(path) for path in discovered}
    if len(faa_by_mag) != len(discovered):
        raise ValueError("Duplicate .faa basenames are not supported in Bakta input.")

bins = Path(snakemake.output.bins)
if bins.exists():
    shutil.rmtree(bins)
bins.mkdir(parents=True, exist_ok=True)
Path(snakemake.log[0]).parent.mkdir(parents=True, exist_ok=True)
extension = str(snakemake.params.extension).lstrip(".")
with open(snakemake.log[0], "w", encoding="utf-8") as log_handle:
    for mag, faa in sorted(faa_by_mag.items()):
        staged = bins / f"{mag}.{extension}"
        staged.symlink_to(Path(faa).resolve())
        log_handle.write(f"{mag}\t{faa}\t{staged}\n")
