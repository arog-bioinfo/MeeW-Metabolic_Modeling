#!/usr/bin/env python3
"""Download and decompress version-pinned integration FASTA fixtures."""

import argparse
import csv
import gzip
import hashlib
import os
import shutil
import tempfile
import urllib.request
from pathlib import Path

MANIFEST_COLUMNS = ["fixture_id", "accession", "source_url", "archive_md5"]
COPY_CHUNK_SIZE = 1024 * 1024


def file_digest(path, algorithm="md5"):
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(COPY_CHUNK_SIZE), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_manifest(path):
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != MANIFEST_COLUMNS:
            raise ValueError(
                f"Manifest columns must be exactly {MANIFEST_COLUMNS}; "
                f"found {reader.fieldnames}."
            )
        records = list(reader)

    if not records:
        raise ValueError("Fixture manifest is empty.")

    seen_ids = set()
    for line_number, record in enumerate(records, start=2):
        if any(not record[column].strip() for column in MANIFEST_COLUMNS):
            raise ValueError(f"Manifest line {line_number} has an empty field.")
        fixture_id = record["fixture_id"]
        if fixture_id in seen_ids:
            raise ValueError(f"Duplicate fixture ID in manifest: {fixture_id}")
        if Path(fixture_id).name != fixture_id or fixture_id in {".", ".."}:
            raise ValueError(f"Unsafe fixture ID in manifest: {fixture_id}")
        expected_md5 = record["archive_md5"].lower()
        if len(expected_md5) != 32 or any(
            character not in "0123456789abcdef" for character in expected_md5
        ):
            raise ValueError(
                f"Invalid archive MD5 on manifest line {line_number}: {expected_md5}"
            )
        record["archive_md5"] = expected_md5
        seen_ids.add(fixture_id)
    return records


def verify_archive(path, expected_md5):
    observed_md5 = file_digest(path)
    if observed_md5 != expected_md5:
        raise ValueError(
            f"Checksum mismatch for {path}: expected {expected_md5}, "
            f"observed {observed_md5}. Remove the corrupt file before retrying."
        )
    return observed_md5


def download_atomic(url, destination, expected_md5):
    request = urllib.request.Request(
        url, headers={"User-Agent": "MeeW-Metabolic-Modeling-integration-test/1"}
    )
    temporary_path = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="wb",
            prefix=f".{destination.name}.",
            suffix=".part",
            dir=destination.parent,
            delete=False,
        ) as temporary:
            temporary_path = Path(temporary.name)
            with urllib.request.urlopen(request, timeout=120) as response:
                shutil.copyfileobj(response, temporary, COPY_CHUNK_SIZE)
            temporary.flush()
            os.fsync(temporary.fileno())

        verify_archive(temporary_path, expected_md5)
        os.replace(temporary_path, destination)
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)


def decompressed_digest(archive):
    digest = hashlib.sha256()
    with gzip.open(archive, "rb") as source:
        for chunk in iter(lambda: source.read(COPY_CHUNK_SIZE), b""):
            digest.update(chunk)
    return digest.hexdigest()


def decompress_atomic(archive, destination):
    if destination.exists():
        expected_digest = decompressed_digest(archive)
        observed_digest = file_digest(destination, algorithm="sha256")
        if observed_digest != expected_digest:
            raise ValueError(
                f"Existing FASTA does not match the verified archive: {destination}. "
                "Remove the corrupt file before retrying."
            )
        return "verified existing"

    temporary_path = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="wb",
            prefix=f".{destination.name}.",
            suffix=".part",
            dir=destination.parent,
            delete=False,
        ) as temporary:
            temporary_path = Path(temporary.name)
            with gzip.open(archive, "rb") as source:
                shutil.copyfileobj(source, temporary, COPY_CHUNK_SIZE)
            temporary.flush()
            os.fsync(temporary.fileno())
        os.replace(temporary_path, destination)
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)
    return "decompressed"


def prepare_fixture(record, output_dir):
    fixture_id = record["fixture_id"]
    archive = output_dir / f"{fixture_id}.faa.gz"
    fasta = output_dir / f"{fixture_id}.faa"

    if archive.exists():
        archive_status = "verified existing"
        observed_md5 = verify_archive(archive, record["archive_md5"])
    else:
        download_atomic(record["source_url"], archive, record["archive_md5"])
        observed_md5 = verify_archive(archive, record["archive_md5"])
        archive_status = "downloaded and verified"

    # The compressed payload is verified above before any decompression occurs.
    fasta_status = decompress_atomic(archive, fasta)
    print(
        f"{fixture_id}: archive {archive_status}; MD5 {observed_md5}; "
        f"FASTA {fasta_status}: {fasta}"
    )


def parse_args():
    script_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        type=Path,
        default=script_dir / "fixture_manifest.tsv",
        help="Tracked fixture manifest (default: beside this script).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=script_dir / "fixtures",
        help="Ignored destination for downloaded archives and FASTAs.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    records = load_manifest(args.manifest.resolve())
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    for record in records:
        prepare_fixture(record, output_dir)


if __name__ == "__main__":
    main()
