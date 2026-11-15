"""Focused input-selection and scratch-staging helpers for CarveMe."""

import csv
import os
import shutil
import tempfile
from contextlib import contextmanager
from pathlib import Path
from typing import NamedTuple

VALID_UNIVERSES = frozenset({"bacteria", "archaea"})


class UniverseSelection(NamedTuple):
    """One resolved CarveMe universe and its concise provenance."""

    universe: str
    source: str
    source_path: str
    reason: str


def _expanded_path(value):
    """Expand environment variables and ``~`` without resolving relative paths."""
    return os.path.expanduser(os.path.expandvars(str(value)))


def _checked_universe(value, *, mag, line_number):
    universe = "" if value is None else str(value).strip()
    if not universe:
        return None
    if universe not in VALID_UNIVERSES:
        raise ValueError(
            "FAA sample TSV has invalid universe "
            f"{universe!r} for MAG {mag!r} on line {line_number}; expected "
            "'bacteria', 'archaea', or a blank cell."
        )
    return universe


def read_faa_samples(path, expand_path=_expanded_path):
    """Read the small standalone FAA sheet with its optional universe column."""
    records = {}
    path = Path(path)
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t", strict=True)
        fieldnames = reader.fieldnames or []
        if any(not field for field in fieldnames) or len(fieldnames) != len(
            set(fieldnames)
        ):
            raise ValueError(
                f"FAA sample TSV has invalid or duplicate columns: {path}."
            )
        missing = [column for column in ("mag", "path") if column not in fieldnames]
        if missing:
            raise ValueError(
                "FAA sample TSV is missing required column(s): "
                f"{', '.join(missing)}."
            )
        has_universe = "universe" in fieldnames
        for line_number, row in enumerate(reader, start=2):
            if None in row or any(value is None for value in row.values()):
                raise ValueError(
                    f"FAA sample TSV has a malformed row at {path}:{line_number}."
                )
            mag = row["mag"].strip()
            faa_path = row["path"].strip()
            if not mag or not faa_path:
                raise ValueError(
                    "FAA sample TSV contains empty mag/path value on line "
                    f"{line_number}."
                )
            if mag in records:
                raise ValueError(
                    f"FAA sample TSV contains duplicate mag {mag!r} on line "
                    f"{line_number}."
                )
            universe = _checked_universe(
                row.get("universe", ""), mag=mag, line_number=line_number
            )
            records[mag] = {
                "path": expand_path(faa_path),
                "universe": universe,
                "universe_fallback_reason": (
                    "blank universe cell in FAA sample TSV"
                    if has_universe and universe is None
                    else (
                        "FAA sample TSV has no universe column"
                        if not has_universe
                        else ""
                    )
                ),
            }
    return records


def read_gtdb_domains(path):
    """Read Annotation's consolidated GTDB-Tk domain table strictly."""
    path = Path(path)
    try:
        handle = path.open(newline="", encoding="utf-8")
    except OSError as error:
        raise FileNotFoundError(
            f"Requested GTDB-Tk classification TSV is missing or unreadable: "
            f"{path}: {error}"
        ) from error

    domains = {}
    lines = {}
    with handle:
        reader = csv.DictReader(handle, delimiter="\t", strict=True)
        fieldnames = reader.fieldnames or []
        if any(not field for field in fieldnames) or len(fieldnames) != len(
            set(fieldnames)
        ):
            raise ValueError(
                "GTDB-Tk classification TSV has invalid or duplicate columns: "
                f"{path}."
            )
        missing = [
            column for column in ("user_genome", "domain") if column not in fieldnames
        ]
        if missing:
            raise ValueError(
                "GTDB-Tk classification TSV is missing required column(s) "
                f"{', '.join(missing)}: {path}."
            )
        for line_number, row in enumerate(reader, start=2):
            if None in row or any(value is None for value in row.values()):
                raise ValueError(
                    "GTDB-Tk classification TSV has a malformed row at "
                    f"{path}:{line_number}."
                )
            genome_id = row["user_genome"].strip()
            domain = row["domain"].strip()
            if not genome_id:
                raise ValueError(
                    "GTDB-Tk classification TSV has a blank user_genome at "
                    f"{path}:{line_number}."
                )
            if genome_id in domains:
                kind = "conflicting" if domain != domains[genome_id] else "duplicate"
                raise ValueError(
                    f"GTDB-Tk classification TSV contains {kind} rows for "
                    f"user_genome {genome_id!r} at lines {lines[genome_id]} and "
                    f"{line_number}: {path}."
                )
            domains[genome_id] = domain
            lines[genome_id] = line_number
    return domains


def resolve_universe(
    mag,
    *,
    explicit_universe=None,
    explicit_source="",
    explicit_fallback_reason="",
    gtdb_classification="",
):
    """Resolve explicit choice, then exact GTDB domain, then bacteria fallback."""
    mag = str(mag).strip()
    if not mag:
        raise ValueError("Cannot resolve a CarveMe universe for a blank MAG ID.")

    if explicit_universe is not None and str(explicit_universe).strip():
        universe = _checked_universe(explicit_universe, mag=mag, line_number="unknown")
        return UniverseSelection(universe, "explicit TSV", str(explicit_source), "")

    if gtdb_classification:
        gtdb_path = str(gtdb_classification)
        domains = read_gtdb_domains(gtdb_path)
        if mag not in domains:
            return UniverseSelection(
                "bacteria",
                "GTDB fallback",
                gtdb_path,
                f"no exact GTDB-Tk user_genome row for MAG {mag!r}",
            )
        domain = domains[mag]
        if domain in VALID_UNIVERSES:
            return UniverseSelection(domain, "GTDB", gtdb_path, "")
        reason = (
            "blank GTDB-Tk domain"
            if not domain
            else f"unsupported GTDB-Tk domain {domain!r}"
        )
        return UniverseSelection("bacteria", "GTDB fallback", gtdb_path, reason)

    return UniverseSelection(
        "bacteria",
        "default",
        str(explicit_source) if explicit_fallback_reason else "",
        explicit_fallback_reason
        or "no explicit universe or GTDB-Tk classification input",
    )


@contextmanager
def staged_faa(source_faa, tmpdir):
    """Copy an FAA into unique scratch and remove that scratch afterward."""
    source_path = Path(source_faa)
    if not source_path.is_file():
        raise FileNotFoundError(
            f"CarveMe source FAA is not a regular file: {source_path}"
        )

    source = source_path.resolve()
    scratch_root = Path(tmpdir).resolve()
    scratch_root.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="carveme-", dir=scratch_root) as workdir:
        staged = Path(workdir) / source_path.name
        shutil.copyfile(source, staged)
        yield staged
