import json
import stat
import subprocess
import sys
from pathlib import Path

import pytest
import yaml
from jsonschema import Draft202012Validator, ValidationError

SCRIPTS = Path(__file__).resolve().parents[1] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

from common.smetana import (  # noqa: E402
    build_smetana_argv,
    normalize_smetana_config,
    run_smetana_mode,
    selected_smetana_modes,
    validate_media_configuration,
    validate_media_database,
)


def smetana_config(global_mode=True, detailed_mode=False, tool_enabled=True):
    return {
        "input": {},
        "solver": "scip",
        "community": {
            "enabled": True,
            "tools": {"smetana": tool_enabled},
        },
        "smetana": {
            "modes": {
                "global": global_mode,
                "detailed": detailed_mode,
            },
            "media_db": "",
            "media": [],
        },
    }


def load_schema():
    schema_path = (
        Path(__file__).resolve().parents[1]
        / "workflow"
        / "schemas"
        / "config.schema.yaml"
    )
    return yaml.safe_load(schema_path.read_text(encoding="utf-8"))


@pytest.mark.parametrize(
    ("global_mode", "detailed_mode", "tool_enabled", "expected"),
    [
        (True, False, True, ("global",)),
        (False, True, True, ("detailed",)),
        (True, True, True, ("global", "detailed")),
    ],
)
def test_mode_selection_constructs_only_enabled_targets(
    global_mode, detailed_mode, tool_enabled, expected
):
    config = smetana_config(global_mode, detailed_mode, tool_enabled)
    assert selected_smetana_modes(config) == expected


def test_default_mode_is_global_only():
    config = {"community": {"enabled": True, "tools": {"smetana": True}}}
    assert normalize_smetana_config(config) == {
        "modes": {"global": True, "detailed": False},
        "media_db": "",
        "media": [],
    }
    assert selected_smetana_modes(config) == ("global",)


def test_enabled_smetana_rejects_no_mode_but_disabled_tool_does_not():
    config = smetana_config(False, False, True)
    with pytest.raises(ValueError, match="neither smetana.modes.global"):
        selected_smetana_modes(config)

    config["community"]["tools"]["smetana"] = False
    assert selected_smetana_modes(config) == ()


def test_schema_enforces_modes_and_media_pairing():
    validator = Draft202012Validator(load_schema())
    validator.validate(smetana_config(True, False, True))

    no_mode = smetana_config(False, False, True)
    with pytest.raises(ValidationError):
        validator.validate(no_mode)

    no_mode["community"]["tools"]["smetana"] = False
    validator.validate(no_mode)

    database_only = smetana_config()
    database_only["smetana"]["media_db"] = "media.tsv"
    with pytest.raises(ValidationError):
        validator.validate(database_only)

    media_only = smetana_config()
    media_only["smetana"]["media"] = ["M9"]
    with pytest.raises(ValidationError):
        validator.validate(media_only)

    paired = smetana_config()
    paired["smetana"].update(media_db="media.tsv", media=["M9"])
    validator.validate(paired)


def test_media_database_validation(tmp_path):
    media_db = tmp_path / "media database.tsv"
    media_db.write_text("medium\tcompound\nM9\tglc__D\nLB\tala__L\n", encoding="utf-8")
    assert validate_media_database(media_db, ["M9", "LB"]) == ("M9", "LB")

    with pytest.raises(ValueError, match="absent"):
        validate_media_database(media_db, ["unknown"])

    missing_column = tmp_path / "missing.tsv"
    missing_column.write_text("name\tcompound\nM9\tglc__D\n", encoding="utf-8")
    with pytest.raises(ValueError, match="'medium' column"):
        validate_media_database(missing_column, ["M9"])


def test_media_db_and_list_must_be_paired():
    assert validate_media_configuration("", []) == ("", ())
    assert validate_media_configuration("media.tsv", ["M9"]) == (
        "media.tsv",
        ("M9",),
    )
    with pytest.raises(ValueError, match="must be supplied together"):
        validate_media_configuration("media.tsv", [])
    with pytest.raises(ValueError, match="must be supplied together"):
        validate_media_configuration("", ["M9"])


def test_argument_construction_for_default_and_configured_media():
    default = build_smetana_argv(
        "/tool/smetana",
        ["model one.xml", "model-two.xml"],
        "global",
        "scip",
        "smetana",
    )
    assert default == [
        "/tool/smetana",
        "model one.xml",
        "model-two.xml",
        "-g",
        "--solver",
        "scip",
        "--flavor",
        "fbc2",
        "-v",
        "-o",
        "smetana",
    ]
    assert "--mediadb" not in default
    assert "-m" not in default

    configured = build_smetana_argv(
        "/tool/smetana",
        ["model.xml"],
        "detailed",
        "cplex",
        "smetana",
        "media db.tsv",
        ["M9", "LB"],
    )
    assert "-d" in configured
    assert configured[configured.index("--mediadb") + 1] == "media db.tsv"
    assert configured[configured.index("-m") + 1] == "M9,LB"
    assert configured.count("-m") == 1

    with pytest.raises(ValueError, match="Unsupported SMETANA mode"):
        build_smetana_argv("smetana", ["model.xml"], "combined", "scip", "x")


def write_fake_smetana(path, fail_detailed=False):
    path.write_text(
        """#!/usr/bin/env python3
import json
import os
import pathlib
import sys

args = sys.argv[1:]
mode = "global" if "-g" in args else "detailed"
prefix = args[args.index("-o") + 1]
mode_index = args.index("-g") if "-g" in args else args.index("-d")
models = args[:mode_index]
payload = {
    "cwd": os.getcwd(),
    "args": args,
    "models": {name: pathlib.Path(name).read_text() for name in models},
}
if "--mediadb" in args:
    media_db = args[args.index("--mediadb") + 1]
    payload["media_db"] = pathlib.Path(media_db).read_text()
pathlib.Path(f"{prefix}_{mode}.tsv").write_text(json.dumps(payload) + "\\n")
pathlib.Path("unrelated.tsv").write_text("must not be published\\n")
print(f"fake SMETANA {mode}")
if mode == "detailed" and FAIL_DETAILED:
    raise SystemExit(23)
""".replace("FAIL_DETAILED", repr(fail_detailed)),
        encoding="utf-8",
    )
    path.chmod(path.stat().st_mode | stat.S_IXUSR)


def test_scratch_staging_cleanup_and_exact_atomic_publication(tmp_path):
    inputs = tmp_path / "source inputs"
    inputs.mkdir()
    model_one = inputs / "model one.xml"
    model_two = inputs / "model-two.xml"
    model_one.write_text("model one source\n", encoding="utf-8")
    model_two.write_text("model two source\n", encoding="utf-8")
    media_db = inputs / "custom media.tsv"
    media_db.write_text("medium\tcompound\nM9\tglc__D\n", encoding="utf-8")
    source_contents = {
        path: path.read_bytes() for path in [model_one, model_two, media_db]
    }

    executable = tmp_path / "fake smetana"
    write_fake_smetana(executable)
    scratch_root = tmp_path / "scratch root"
    output = tmp_path / "published results" / "global.tsv"
    log = tmp_path / "logs" / "global.log"
    log.parent.mkdir()
    log.write_text("CPLEX preflight preserved\n", encoding="utf-8")

    run_smetana_mode(
        executable=executable,
        models=[model_one, model_two],
        mode="global",
        solver="scip",
        output=output,
        log=log,
        tmpdir=scratch_root,
        media_db=media_db,
        media=["M9"],
    )

    payload = json.loads(output.read_text(encoding="utf-8"))
    scratch_job = Path(payload["cwd"])
    assert scratch_job.parent == scratch_root
    assert not scratch_job.exists()
    assert list(scratch_root.iterdir()) == []
    assert payload["models"] == {
        "model one.xml": "model one source\n",
        "model-two.xml": "model two source\n",
    }
    assert payload["media_db"] == "medium\tcompound\nM9\tglc__D\n"
    assert payload["args"][payload["args"].index("-m") + 1] == "M9"
    assert not (output.parent / "unrelated.tsv").exists()
    assert list(output.parent.glob(f".{output.name}.*.tmp")) == []
    assert log.read_text(encoding="utf-8").splitlines() == [
        "CPLEX preflight preserved",
        "fake SMETANA global",
    ]
    assert {path: path.read_bytes() for path in source_contents} == source_contents


def test_failed_detailed_mode_cannot_remove_completed_global_output(tmp_path):
    model = tmp_path / "model.xml"
    model.write_text("unchanged model\n", encoding="utf-8")
    executable = tmp_path / "fake-smetana"
    write_fake_smetana(executable, fail_detailed=True)
    scratch_root = tmp_path / "scratch"
    global_output = tmp_path / "results" / "global.tsv"
    detailed_output = tmp_path / "results" / "detailed.tsv"

    run_smetana_mode(
        executable=executable,
        models=[model],
        mode="global",
        solver="scip",
        output=global_output,
        log=tmp_path / "global.log",
        tmpdir=scratch_root,
    )
    completed_global = global_output.read_bytes()

    with pytest.raises(
        subprocess.CalledProcessError, match="returned non-zero exit status 23"
    ):
        run_smetana_mode(
            executable=executable,
            models=[model],
            mode="detailed",
            solver="scip",
            output=detailed_output,
            log=tmp_path / "detailed.log",
            tmpdir=scratch_root,
        )

    assert global_output.read_bytes() == completed_global
    assert not detailed_output.exists()
    assert list(scratch_root.iterdir()) == []
    assert model.read_text(encoding="utf-8") == "unchanged model\n"


def test_carveme_real_command_requests_fbc2():
    carveme_script = SCRIPTS / "carveme.py"
    assert '"--fbc2 "' in carveme_script.read_text(encoding="utf-8")
