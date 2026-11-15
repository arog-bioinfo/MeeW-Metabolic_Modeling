"""Helpers for inspecting collections of workflow paths."""

from collections import Counter
from pathlib import Path


def duplicate_path_stems(paths):
    """Return sorted duplicate path stems in linear rather than quadratic time."""
    stem_counts = Counter(Path(path).stem for path in paths)
    return sorted(stem for stem, count in stem_counts.items() if count > 1)
