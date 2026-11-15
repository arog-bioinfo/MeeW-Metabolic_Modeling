"""Resolve exact Annotation-owned Bakta directory inputs at job time."""

from pathlib import Path

ANNOTATION_FAA_NAME = "annotation.faa"


def annotation_faa(directory):
    """Return the fixed Bakta FAA child after validating its declared parent."""
    parent = Path(directory)
    if not parent.is_dir():
        raise FileNotFoundError(
            f"Declared Annotation Bakta directory does not exist: {parent}"
        )
    faa = parent / ANNOTATION_FAA_NAME
    if not faa.is_file():
        raise FileNotFoundError(
            f"Declared Annotation Bakta directory '{parent}' is missing its fixed "
            f"'{ANNOTATION_FAA_NAME}' child."
        )
    return str(faa)


def annotation_faa_by_mag(mags, directories):
    """Pair caller-supplied MAG IDs with declared Bakta directories in order."""
    mag_ids = [str(mag) for mag in mags]
    parents = [str(directory) for directory in directories]
    if len(mag_ids) != len(parents):
        raise ValueError(
            "Direct Annotation input mismatch: received "
            f"{len(mag_ids)} MAG IDs and {len(parents)} Bakta directories."
        )
    if len(set(mag_ids)) != len(mag_ids):
        raise ValueError("Direct Annotation input contains duplicate MAG IDs.")
    return {
        mag: annotation_faa(parent)
        for mag, parent in zip(mag_ids, parents, strict=True)
    }
