"""Failure-safe lifecycle for CarveMe's bundled BiGG DIAMOND index."""

from __future__ import annotations

import configparser
import fcntl
import hashlib
import importlib.metadata
import json
import os
import shlex
import shutil
import stat
import subprocess
import tempfile
import uuid
from contextlib import contextmanager
from pathlib import Path
from typing import NamedTuple


class CarveMeIndexError(RuntimeError):
    """CarveMe's configured bundled index is missing, invalid, or changed."""


class BundledIndex(NamedTuple):
    """Native CarveMe paths resolved from its installed package configuration."""

    package: Path
    config: Path
    fasta: Path
    database: Path


_WORKFLOW_INVOCATION_ATTRIBUTE = "_meew_carveme_index_invocation"


def workflow_invocation_token(workflow_state, *, token_factory=None):
    """Return one token per Workflow object, stable across module reevaluation."""
    try:
        return getattr(workflow_state, _WORKFLOW_INVOCATION_ATTRIBUTE)
    except AttributeError:
        token = (token_factory or (lambda: uuid.uuid4().hex))()
        setattr(workflow_state, _WORKFLOW_INVOCATION_ATTRIBUTE, token)
        return token


def _resolved(path):
    return Path(path).expanduser().resolve()


def _append_log(log_path, message):
    if not log_path:
        return
    path = _resolved(log_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(f"[MeeW CarveMe index] {message}\n")


def discover_bundled_index(distribution_name="carveme"):
    """Resolve CarveMe's real bundled FASTA and index without importing it."""
    try:
        distribution = importlib.metadata.distribution(distribution_name)
    except importlib.metadata.PackageNotFoundError as error:
        raise CarveMeIndexError(
            "Cannot locate the installed CarveMe distribution in the active "
            "modeling environment."
        ) from error

    package = _resolved(distribution.locate_file("carveme"))
    config_path = package / "config.cfg"
    if not config_path.is_file():
        raise CarveMeIndexError(
            f"CarveMe's installed native configuration is missing: {config_path}"
        )

    parser = configparser.ConfigParser()
    try:
        with config_path.open(encoding="utf-8") as handle:
            parser.read_file(handle)
        fasta_setting = parser.get("generated", "fasta_file")
        database_setting = parser.get("generated", "diamond_db")
    except (OSError, configparser.Error) as error:
        raise CarveMeIndexError(
            f"Cannot read CarveMe's installed native configuration {config_path}: "
            f"{error}"
        ) from error

    fasta = _resolved(package / fasta_setting)
    database = _resolved(package / database_setting)
    if database.suffix != ".dmnd":
        raise CarveMeIndexError(
            "CarveMe's configured DIAMOND database does not end in .dmnd: "
            f"{database}"
        )
    return BundledIndex(package, config_path, fasta, database)


def _regular_nonempty_identity(path, label):
    path = Path(path)
    try:
        details = path.stat()
    except OSError as error:
        raise CarveMeIndexError(
            f"{label} is missing or unreadable: {path}: {error}"
        ) from error
    if not stat.S_ISREG(details.st_mode) or details.st_size <= 0:
        raise CarveMeIndexError(f"{label} is not a non-empty regular file: {path}")
    if not os.access(path, os.R_OK):
        raise CarveMeIndexError(f"{label} is not readable: {path}")
    return {
        "path": str(_resolved(path)),
        "device": details.st_dev,
        "inode": details.st_ino,
        "size": details.st_size,
        "mtime_ns": details.st_mtime_ns,
    }


def _run_logged(argv, log_path):
    argv = [str(argument) for argument in argv]
    _append_log(log_path, f"$ {shlex.join(argv)}")
    try:
        if log_path:
            path = _resolved(log_path)
            path.parent.mkdir(parents=True, exist_ok=True)
            with path.open("a", encoding="utf-8") as handle:
                return subprocess.run(
                    argv,
                    check=False,
                    stdout=handle,
                    stderr=subprocess.STDOUT,
                    text=True,
                )
        return subprocess.run(
            argv,
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
    except OSError as error:
        raise CarveMeIndexError(
            f"Unable to execute native DIAMOND command {shlex.join(argv)}: {error}"
        ) from error


def validate_index(database, *, diamond="diamond", log_path=None):
    """Validate one complete index with native DIAMOND and return its identity."""
    identity = _regular_nonempty_identity(database, "CarveMe DIAMOND index")
    result = _run_logged([diamond, "dbinfo", "--db", database], log_path)
    if result.returncode != 0:
        detail = f" Output: {(result.stdout or '').strip()}" if not log_path else ""
        raise CarveMeIndexError(
            f"Native DIAMOND validation failed for {database} with exit code "
            f"{result.returncode}.{detail}"
        )
    return identity


def cache_lock_path(database):
    """Return the canonical cache-keyed lock shared by every MeeW workflow."""
    database = _resolved(database)
    key = hashlib.sha256(str(database).encode()).hexdigest()[:20]
    return database.parent / f".meew-carveme-index-{key}.lock"


@contextmanager
def _exclusive_cache_lock(database):
    lock_path = cache_lock_path(database)
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with lock_path.open("a+", encoding="utf-8") as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        try:
            yield lock_path
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


@contextmanager
def _private_staging_directory(database):
    database = _resolved(database)
    database.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(
        tempfile.mkdtemp(prefix=".meew-carveme-index-stage-", dir=database.parent)
    )
    try:
        yield staging
    finally:
        if staging.exists():
            shutil.rmtree(staging)


def _path_lexists(path):
    return os.path.lexists(os.fspath(path))


def ensure_bundled_index(paths, *, diamond="diamond", log_path=None):
    """Reuse or atomically build CarveMe's configured single-file index."""
    fasta_identity = _regular_nonempty_identity(
        paths.fasta, "CarveMe bundled BiGG protein FASTA"
    )
    config_identity = _regular_nonempty_identity(
        paths.config, "CarveMe native configuration"
    )

    try:
        database_identity = validate_index(
            paths.database, diamond=diamond, log_path=log_path
        )
    except CarveMeIndexError:
        database_identity = None
    else:
        _append_log(log_path, f"Reused valid bundled index read-only: {paths.database}")
        return "reused", config_identity, fasta_identity, database_identity

    with _exclusive_cache_lock(paths.database):
        validation_error = None
        try:
            database_identity = validate_index(
                paths.database, diamond=diamond, log_path=log_path
            )
        except CarveMeIndexError as error:
            validation_error = error
            database_identity = None
        else:
            _append_log(
                log_path,
                f"Bundled index became ready while waiting; reused read-only: "
                f"{paths.database}",
            )
            return "reused", config_identity, fasta_identity, database_identity

        if _path_lexists(paths.database):
            message = (
                "CarveMe's configured DIAMOND index exists but is invalid. "
                "The pre-existing index was preserved unchanged and will not be "
                f"silently replaced: {paths.database}. Validation error: "
                f"{validation_error}"
            )
            _append_log(log_path, message)
            raise CarveMeIndexError(message) from validation_error

        with _private_staging_directory(paths.database) as staging:
            staged_prefix = staging / paths.database.name.removesuffix(".dmnd")
            staged_database = Path(f"{staged_prefix}.dmnd")
            result = _run_logged(
                [diamond, "makedb", "--in", paths.fasta, "-d", staged_prefix],
                log_path,
            )
            if result.returncode != 0:
                detail = (
                    f" Output: {(result.stdout or '').strip()}" if not log_path else ""
                )
                raise CarveMeIndexError(
                    "DIAMOND failed while building CarveMe's private staged index "
                    f"with exit code {result.returncode}.{detail}"
                )
            validate_index(staged_database, diamond=diamond, log_path=log_path)
            if _path_lexists(paths.database):
                raise CarveMeIndexError(
                    "Refusing to replace a CarveMe index that appeared outside the "
                    f"shared MeeW lock: {paths.database}"
                )
            with staged_database.open("rb") as handle:
                os.fsync(handle.fileno())
            os.replace(staged_database, paths.database)

        database_identity = _regular_nonempty_identity(
            paths.database, "CarveMe DIAMOND index"
        )
        _append_log(
            log_path,
            f"Built, validated, and atomically published bundled index: "
            f"{paths.database}",
        )
        return "prepared", config_identity, fasta_identity, database_identity


def write_readiness(
    path,
    *,
    paths,
    status,
    config_identity,
    fasta_identity,
    database_identity,
    initializer_log,
):
    """Atomically record readiness only after native validation succeeds."""
    destination = _resolved(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{destination.name}.", suffix=".tmp", dir=destination.parent
    )
    temporary = Path(temporary_name)
    payload = {
        "schema_version": 1,
        "status": status,
        "package": str(paths.package),
        "initializer_log": str(initializer_log) if initializer_log else "",
        "config": config_identity,
        "fasta": fasta_identity,
        "database": database_identity,
    }
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, destination)
    finally:
        if temporary.exists():
            temporary.unlink()


def verify_readiness(path, *, paths=None):
    """Reject stale readiness if the active environment/cache identity changed."""
    readiness = _resolved(path)
    try:
        payload = json.loads(readiness.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise CarveMeIndexError(
            f"CarveMe index readiness is missing or malformed: {readiness}: {error}"
        ) from error
    if payload.get("schema_version") != 1 or payload.get("status") not in {
        "prepared",
        "reused",
    }:
        raise CarveMeIndexError(
            f"CarveMe index readiness has an unsupported state: {readiness}"
        )

    paths = paths or discover_bundled_index()
    expected_paths = {
        "config": paths.config,
        "fasta": paths.fasta,
        "database": paths.database,
    }
    for name, expected_path in expected_paths.items():
        recorded = payload.get(name)
        if not isinstance(recorded, dict):
            raise CarveMeIndexError(
                f"CarveMe index readiness is missing {name} identity: {readiness}"
            )
        current = _regular_nonempty_identity(expected_path, f"CarveMe {name}")
        if recorded != current:
            raise CarveMeIndexError(
                f"CarveMe {name} changed after initialization; refusing to launch "
                f"native CarveMe. Re-run the workflow initializer. Recorded "
                f"{recorded!r}; current {current!r}. Initialization log: "
                f"{payload.get('initializer_log', '') or '(not recorded)'}"
            )
    return payload
