from copy import deepcopy
from pathlib import Path

import yaml
from jsonschema import Draft202012Validator

ROOT = Path(__file__).resolve().parents[1]
WORKFLOW = ROOT / ".github" / "workflows" / "main.yaml"
DISABLED_WORKFLOW = ROOT / ".github" / "workflows-disabled" / "main.yaml"
BASE_CONFIG = ROOT / ".test" / "config" / "config.yaml"
CHECKM2_CONFIG = ROOT / ".test" / "config" / "config.automatic-checkm2.yaml"
CONFIG_SCHEMA = ROOT / "workflow" / "schemas" / "config.schema.yaml"

EXPECTED_SUITES = [
    {"name": "external QA fixture", "config": ".test/config/config.yaml"},
    {
        "name": "automatic CheckM2 first QA",
        "config": ".test/config/config.automatic-checkm2.yaml",
    },
]


def _load_yaml(path, loader=yaml.BaseLoader):
    return yaml.load(path.read_text(encoding="utf-8"), Loader=loader)


def _named_step(job, name):
    return next(step for step in job["steps"] if step.get("name") == name)


def _merge_config(base, override):
    merged = deepcopy(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = _merge_config(merged[key], value)
        else:
            merged[key] = deepcopy(value)
    return merged


def test_template_plus_ci_contract():
    workflow = _load_yaml(WORKFLOW)

    assert WORKFLOW.is_file()
    assert not DISABLED_WORKFLOW.exists()
    triggers = workflow["on"]
    assert set(triggers) == {"push", "pull_request", "workflow_dispatch"}
    assert triggers["push"]["branches"] == ["main"]
    assert triggers["pull_request"]["branches"] == ["main"]
    assert workflow["permissions"] == {"contents": "read"}
    assert workflow["concurrency"] == {
        "group": "${{ github.workflow }}-${{ github.ref }}",
        "cancel-in-progress": "true",
    }

    jobs = workflow["jobs"]
    assert set(jobs) == {"formatting", "linting", "python-tests", "testing"}
    assert "needs" not in jobs["formatting"]
    assert "needs" not in jobs["linting"]
    assert "needs" not in jobs["python-tests"]
    assert jobs["testing"]["needs"] == [
        "formatting",
        "linting",
        "python-tests",
    ]
    assert {name: job["timeout-minutes"] for name, job in jobs.items()} == {
        "formatting": "15",
        "linting": "15",
        "python-tests": "20",
        "testing": "20",
    }
    for job in jobs.values():
        checkout_refs = [
            step["uses"]
            for step in job["steps"]
            if step.get("uses", "").startswith("actions/checkout@")
        ]
        assert checkout_refs == ["actions/checkout@v6"]

    formatting = jobs["formatting"]
    assert formatting["steps"][0] == {
        "uses": "actions/checkout@v6",
        "with": {"fetch-depth": "0"},
    }
    formatter = _named_step(formatting, "Check formatting")
    assert formatter["uses"] == "super-linter/super-linter@v8"
    assert {
        key: value
        for key, value in formatter["env"].items()
        if key.startswith("VALIDATE_")
    } == {
        "VALIDATE_ALL_CODEBASE": "true",
        "VALIDATE_PYTHON_BLACK": "true",
        "VALIDATE_SNAKEMAKE_SNAKEFMT": "true",
        "VALIDATE_YAML_PRETTIER": "true",
    }

    lint = _named_step(jobs["linting"], "Lint workflow")
    assert lint["uses"] == "snakemake/snakemake-github-action@v2"
    assert lint["with"] == {
        "directory": ".",
        "snakefile": "workflow/Snakefile",
        "args": "--lint",
        "snakemake-version": "9.19.0",
    }

    python_job = jobs["python-tests"]
    setup_python = next(
        step
        for step in python_job["steps"]
        if step.get("uses") == "actions/setup-python@v6"
    )
    assert setup_python["with"] == {
        "python-version": "3.13",
        "cache": "pip",
        "cache-dependency-path": (
            "tests/requirements.txt\ntests/requirements-workflow.txt\n"
        ),
    }
    assert _named_step(python_job, "Install test dependencies")["run"] == (
        "python -m pip install -r tests/requirements-workflow.txt"
    )
    assert _named_step(python_job, "Run Python tests")["run"] == (
        "python -m pytest -q -ra tests"
    )

    testing = jobs["testing"]
    assert testing["strategy"] == {
        "fail-fast": "false",
        "matrix": {"suite": EXPECTED_SUITES},
    }
    assert testing["steps"][0] == {"uses": "actions/checkout@v6"}
    assert len(testing["steps"]) == 2
    dry_run = _named_step(testing, "Dry-run workflow")
    assert dry_run["uses"] == "snakemake/snakemake-github-action@v2"
    assert dry_run["with"]["directory"] == ".test"
    assert dry_run["with"]["snakefile"] == "workflow/Snakefile"
    assert dry_run["with"]["snakemake-version"] == "9.19.0"
    dag_args = " ".join(dry_run["with"]["args"].split())
    assert dag_args == (
        "--configfile ${{ matrix.suite.config }} "
        "--dry-run --printshellcmds --cores 2"
    )
    for forbidden in (
        "--software-deployment-method",
        "--sdm",
        "conda",
        "create-env",
        "--report",
        "container",
        "artifact",
    ):
        assert forbidden not in str(testing).lower()


def test_ci_dependencies_are_pinned():
    requirements = (ROOT / "tests" / "requirements.txt").read_text(encoding="utf-8")
    workflow_requirements = (ROOT / "tests" / "requirements-workflow.txt").read_text(
        encoding="utf-8"
    )

    assert requirements.splitlines() == [
        "pytest==9.0.2",
        "PyYAML==6.0.3",
        "jsonschema==4.26.0",
        "pandas==2.3.3",
    ]
    assert workflow_requirements.splitlines() == [
        "-r requirements.txt",
        "snakemake==9.19.0",
    ]


def test_automatic_checkm2_override_merges_into_valid_mock_fixture():
    base = _load_yaml(BASE_CONFIG, yaml.SafeLoader)
    override = _load_yaml(CHECKM2_CONFIG, yaml.SafeLoader)
    assert override == {
        "qa": {"tsv": ""},
        "checkm2": {"enabled": True, "download_database": False},
    }

    merged = _merge_config(base, override)
    assert merged["execution"] == {"mock": True}
    assert merged["input"] == {"faa_globs": ["data/mags/*.faa"], "faa_dir": ""}
    assert merged["qa"]["enabled"] is True
    assert merged["qa"]["tsv"] == ""
    assert merged["checkm2"]["enabled"] is True
    assert merged["checkm2"]["download_database"] is False
    assert merged["checkm2"]["database_path"] == "resources/checkm2_db"

    fixture_faas = sorted((ROOT / ".test").glob(merged["input"]["faa_globs"][0]))
    assert [path.name for path in fixture_faas] == ["bin1.faa", "bin2.faa", "bin3.faa"]
    assert all(path.is_file() for path in fixture_faas)

    # common.smk adds these documented defaults before validating the config.
    merged["smetana"] = {
        "modes": {"global": True, "detailed": False},
        "media_db": "",
        "media": [],
    }
    schema = _load_yaml(CONFIG_SCHEMA, yaml.SafeLoader)
    Draft202012Validator(schema).validate(merged)
