"""Shared model-selection helpers for community analysis scripts."""

from pathlib import Path


def model_id_from_path(model_path):
    """Return the model identifier encoded in a model filename."""
    return Path(model_path).stem


def select_memote_passing_models(models, filtered_models_file):
    """Select models listed by the Memote checkpoint, preserving input order."""
    model_paths = [models] if isinstance(models, str) else list(models)
    with Path(filtered_models_file).open(encoding="utf-8") as handle:
        passing_model_ids = {line.strip() for line in handle if line.strip()}
    return [
        model_path
        for model_path in model_paths
        if model_id_from_path(model_path) in passing_model_ids
    ]
