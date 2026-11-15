"""Configuration and execution helpers for one SMETANA mode at a time."""

import csv
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

MODE_FLAGS = {"global": "-g", "detailed": "-d"}


def normalize_smetana_config(config):
    """Populate the documented SMETANA defaults in a workflow config."""
    smetana = config.setdefault("smetana", {})
    modes = smetana.setdefault("modes", {})
    modes.setdefault("global", True)
    modes.setdefault("detailed", False)
    smetana.setdefault("media_db", "")
    smetana.setdefault("media", [])
    return smetana


def selected_smetana_modes(config):
    """Return enabled target modes, respecting the community/tool switches."""
    community = config.get("community", {})
    tool_enabled = bool(community.get("tools", {}).get("smetana", True))
    if not bool(community.get("enabled", True)) or not tool_enabled:
        return ()

    modes = normalize_smetana_config(config)["modes"]
    selected = tuple(mode for mode in MODE_FLAGS if bool(modes.get(mode, False)))
    if not selected:
        raise ValueError(
            "SMETANA is enabled, but neither smetana.modes.global nor "
            "smetana.modes.detailed is true."
        )
    return selected


def validate_media_configuration(media_db, media):
    """Validate and normalize the media DB/list pairing contract."""
    database = str(media_db or "").strip()
    selected = tuple(str(name).strip() for name in (media or []))
    if any(not name for name in selected):
        raise ValueError("smetana.media entries must be nonempty names.")
    if bool(database) != bool(selected):
        raise ValueError(
            "smetana.media_db and a nonempty smetana.media list must be supplied "
            "together, or both must be empty."
        )
    return database, selected


def validate_media_database(media_db, media):
    """Require a medium column and confirm every selected medium is present."""
    database, selected = validate_media_configuration(media_db, media)
    if not database:
        return ()

    with open(database, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if "medium" not in (reader.fieldnames or []):
            raise ValueError(
                f"SMETANA media database {database!r} is missing required "
                "'medium' column."
            )
        available = {
            str(row.get("medium", "")).strip()
            for row in reader
            if str(row.get("medium", "")).strip()
        }

    unknown = sorted(set(selected) - available)
    if unknown:
        raise ValueError(
            f"Selected SMETANA media are absent from {database!r}: {unknown}"
        )
    return selected


def build_smetana_argv(
    executable,
    models,
    mode,
    solver,
    prefix,
    media_db="",
    media=(),
):
    """Build a shell-free SMETANA command for one approved mode."""
    if mode not in MODE_FLAGS:
        raise ValueError(
            f"Unsupported SMETANA mode {mode!r}; expected one of "
            f"{sorted(MODE_FLAGS)}."
        )

    database, selected = validate_media_configuration(media_db, media)
    command = [
        str(executable),
        *(str(model) for model in models),
        MODE_FLAGS[mode],
        "--solver",
        str(solver),
        "--flavor",
        "fbc2",
        "-v",
        "-o",
        str(prefix),
    ]
    if database:
        command.extend(["--mediadb", str(database), "-m", ",".join(selected)])
    return command


def _stage_models(models, workdir):
    staged = []
    used_names = set()
    for model in models:
        source = Path(model).resolve()
        if not source.is_file():
            raise FileNotFoundError(f"SMETANA model does not exist: {source}")
        if source.name in used_names:
            raise ValueError(
                f"SMETANA model basenames must be unique for scratch staging: "
                f"{source.name!r}"
            )
        used_names.add(source.name)
        destination = workdir / source.name
        shutil.copy2(source, destination)
        staged.append(destination.name)
    return staged


def _publish_atomic(source, destination):
    destination = Path(destination).resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = None
    try:
        with tempfile.NamedTemporaryFile(
            prefix=f".{destination.name}.",
            suffix=".tmp",
            dir=destination.parent,
            delete=False,
        ) as temporary:
            temporary_path = Path(temporary.name)
        shutil.copyfile(source, temporary_path)
        os.replace(temporary_path, destination)
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)


def run_smetana_mode(
    *,
    executable,
    models,
    mode,
    solver,
    output,
    log,
    tmpdir,
    media_db="",
    media=(),
):
    """Run one SMETANA mode in unique scratch and atomically publish its output."""
    if mode not in MODE_FLAGS:
        raise ValueError(
            f"Unsupported SMETANA mode {mode!r}; expected one of "
            f"{sorted(MODE_FLAGS)}."
        )
    model_paths = [Path(model).resolve() for model in models]
    if not model_paths:
        raise ValueError("No Memote-passing models are available for SMETANA.")

    database, selected = validate_media_configuration(media_db, media)
    if database:
        database = str(Path(database).resolve())
        validate_media_database(database, selected)

    scratch_root = Path(tmpdir).resolve()
    scratch_root.mkdir(parents=True, exist_ok=True)
    log_path = Path(log).resolve()
    log_path.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(
        prefix=f"smetana-{mode}-", dir=scratch_root
    ) as temporary_directory:
        workdir = Path(temporary_directory)
        staged_models = _stage_models(model_paths, workdir)
        staged_database = ""
        if database:
            staged_database = "media_db.tsv"
            shutil.copy2(database, workdir / staged_database)

        prefix = "smetana"
        command = build_smetana_argv(
            executable,
            staged_models,
            mode,
            solver,
            prefix,
            staged_database,
            selected,
        )
        with log_path.open("a", encoding="utf-8") as log_handle:
            subprocess.run(
                command,
                cwd=workdir,
                stdout=log_handle,
                stderr=subprocess.STDOUT,
                check=True,
            )

        expected = workdir / f"{prefix}_{mode}.tsv"
        if not expected.is_file():
            raise FileNotFoundError(
                f"SMETANA {mode} mode completed without expected output: "
                f"{expected.name}"
            )
        _publish_atomic(expected, output)
