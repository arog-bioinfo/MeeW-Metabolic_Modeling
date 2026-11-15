import csv
import json
import os
import runpy
import shutil
import stat
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
SNAKEFILE = ROOT / "workflow" / "Snakefile"
STAGE_CHECKM2 = ROOT / "workflow" / "scripts" / "stage_checkm2_faa.py"
DIRECT_INPUTS_HELPER = ROOT / "workflow" / "scripts" / "meew_metabolic_direct_inputs.py"


class _NamedValues(dict):
    __getattr__ = dict.__getitem__


def _mock_config(sample, faa=None):
    config = yaml.safe_load((ROOT / "config" / "config.yaml").read_text())
    config["sample"] = sample
    config["input"] = {
        "faa_samples": "",
        "faa_globs": [str(faa)] if faa is not None else [],
        "faa_dir": "",
        "bakta_dir": "",
        "bakta_done": "",
    }
    config["solver"] = "scip"
    config["execution"]["mock"] = True
    config["qa"]["enabled"] = False
    config["community"] = {
        "enabled": False,
        "tools": {"smetana": False, "micom": False, "steadycom": False},
    }
    return config


def _run(workdir, snakefile, configfile, *arguments, environment=None):
    run_environment = os.environ.copy()
    run_environment.update(
        PYTHONDONTWRITEBYTECODE="1",
        PYTHONPYCACHEPREFIX=str(workdir / "pycache"),
        XDG_CACHE_HOME=str(workdir / "xdg-cache"),
    )
    if environment:
        run_environment.update(environment)
    return subprocess.run(
        [
            sys.executable,
            "-m",
            "snakemake",
            "--snakefile",
            str(snakefile),
            "--directory",
            str(workdir),
            "--configfile",
            str(configfile),
            "--nolock",
            "--cores",
            "2",
            *arguments,
        ],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
        env=run_environment,
    )


def _assert_complete_mock_cohort(workdir, sample, mag, results_prefix="results"):
    root = workdir / results_prefix / sample
    assert (root / "qa" / "filtered_mags.txt").read_text(encoding="utf-8") == (
        f"{mag}\n"
    )
    score = root / "memote" / mag / "score.json"
    assert json.loads(score.read_text(encoding="utf-8"))["score"]["total_score"] == 1.0
    with (root / "memote" / "memote_quality.tsv").open(
        newline="", encoding="utf-8"
    ) as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [(row["mag"], row["status"]) for row in rows] == [(mag, "PASS")]
    assert (root / "memote" / "filtered_models.txt").read_text(
        encoding="utf-8"
    ) == f"{mag}\n"


def _write_fake_checkm2(workdir):
    fake_bin = workdir / "fake bin"
    fake_bin.mkdir()
    checkm2 = fake_bin / "checkm2"
    checkm2.write_text(
        f"#!{sys.executable}\n"
        "import sys\n"
        "from pathlib import Path\n"
        "args = sys.argv[1:]\n"
        "assert args[0] == 'predict'\n"
        "input_dir = Path(args[args.index('--input') + 1])\n"
        "output_dir = Path(args[args.index('--output-directory') + 1])\n"
        "output_dir.mkdir(parents=True, exist_ok=True)\n"
        "mags = sorted(path.stem for path in input_dir.glob('*.faa'))\n"
        "(output_dir / 'quality_report.tsv').write_text(\n"
        "    'Name\\tCompleteness\\tContamination\\n'\n"
        "    + ''.join(f'{mag}\\t90\\t1\\n' for mag in mags),\n"
        "    encoding='utf-8',\n"
        ")\n",
        encoding="utf-8",
    )
    checkm2.chmod(checkm2.stat().st_mode | stat.S_IXUSR)
    return fake_bin


def test_standalone_native_checkpoint_helpers_execute_and_cache(tmp_path):
    faa = tmp_path / "standalone-mag.faa"
    faa.write_text(">standalone-mag\nMPEPTIDE\n", encoding="utf-8")
    configfile = tmp_path / "standalone.yaml"
    configfile.write_text(
        yaml.safe_dump(_mock_config("standalone", faa), sort_keys=False),
        encoding="utf-8",
    )

    first = _run(tmp_path, SNAKEFILE, configfile, "all")
    assert first.returncode == 0, first.stdout + first.stderr
    _assert_complete_mock_cohort(tmp_path, "standalone", "standalone-mag")

    cached = _run(tmp_path, SNAKEFILE, configfile, "all")
    assert cached.returncode == 0, cached.stdout + cached.stderr
    assert "Nothing to be done" in cached.stdout + cached.stderr


def test_unbound_standalone_automatic_checkm2_is_unchanged(tmp_path):
    faa = tmp_path / "standalone-auto-mag.faa"
    faa.write_text(">standalone-auto-mag\nMPEPTIDE\n", encoding="utf-8")
    config = _mock_config("standalone-auto", faa)
    config["qa"]["enabled"] = True
    config["qa"]["tsv"] = ""
    config["checkm2"]["database_path"] = str(tmp_path / "mock-checkm2-db")
    configfile = tmp_path / "standalone-auto.yaml"
    configfile.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    fake_bin = _write_fake_checkm2(tmp_path)
    environment = {"PATH": f"{fake_bin}{os.pathsep}{os.environ['PATH']}"}

    result = _run(
        tmp_path,
        SNAKEFILE,
        configfile,
        "--printshellcmds",
        "all",
        environment=environment,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    _assert_complete_mock_cohort(tmp_path, "standalone-auto", "standalone-auto-mag")
    staged = tmp_path / "results" / "standalone-auto" / "qa" / "checkm2" / "input"
    assert (staged / "standalone-auto-mag.faa").is_symlink()
    assert (staged / "standalone-auto-mag.faa").resolve() == faa.resolve()
    assert "meew_metabolic_direct_inputs.py" not in result.stdout + result.stderr


@pytest.mark.parametrize(
    "unsafe_mag",
    [
        "../../escaped-output",
        r"..\escaped-output",
        "nested/escaped-output",
        r"nested\escaped-output",
        "control\nvalue",
        "nul\x00value",
    ],
    ids=[
        "slash-traversal",
        "backslash-traversal",
        "slash-component",
        "backslash-component",
        "control",
        "nul",
    ],
)
def test_bound_unsafe_mag_ids_fail_before_outputs(tmp_path, unsafe_mag):
    configfile = tmp_path / "unsafe.yaml"
    configfile.write_text(
        yaml.safe_dump(_mock_config("unsafe-bound"), sort_keys=False),
        encoding="utf-8",
    )
    harness = tmp_path / "unsafe-bound.smk"
    harness.write_text(
        f"UNSAFE_MAG = {json.dumps(unsafe_mag)}\n\n"
        "module modeling:\n"
        f"    snakefile: {json.dumps(str(SNAKEFILE))}\n"
        "    config: config\n\n"
        "use rule * from modeling exclude all as modeling_*\n\n"
        "modeling.bind_protein_sources(\n"
        "    sources={UNSAFE_MAG: 'future/protein parent'}\n"
        ")\n\n"
        "rule produce_parent:\n"
        "    output:\n"
        "        directory('future/protein parent'),\n"
        "    shell:\n"
        "        \"mkdir -p {output:q}; printf '>unsafe\\nM\\n' > {output:q}/annotation.faa\"\n\n"
        "rule all:\n"
        "    input:\n"
        "        modeling.FINAL_TARGETS,\n"
        "        'results/unsafe-bound/qa/checkm2/input',\n",
        encoding="utf-8",
    )

    result = _run(tmp_path, harness, configfile, "all")
    assert result.returncode != 0
    assert "Invalid bound protein source MAG ID" in result.stdout + result.stderr
    assert not (tmp_path / "results").exists()
    assert not (tmp_path / "future").exists()


def test_source_cached_checkm2_stage_loads_declared_helper(tmp_path):
    cached_script = tmp_path / "source cache" / "scripts" / STAGE_CHECKM2.name
    cached_helper = tmp_path / "source cache" / "assets" / DIRECT_INPUTS_HELPER.name
    cached_script.parent.mkdir(parents=True)
    cached_helper.parent.mkdir(parents=True)
    shutil.copy2(STAGE_CHECKM2, cached_script)
    shutil.copy2(DIRECT_INPUTS_HELPER, cached_helper)

    protein_parent = tmp_path / "future protein parent"
    protein_parent.mkdir()
    faa = protein_parent / "annotation.faa"
    faa.write_text(">cached-mag\nMPEPTIDE\n", encoding="utf-8")
    bins = tmp_path / "staged bins"
    log = tmp_path / "logs" / "stage.log"
    fake_snakemake = SimpleNamespace(
        input=_NamedValues(direct_inputs_helper=str(cached_helper)),
        output=_NamedValues(bins=str(bins)),
        params=_NamedValues(
            faa_by_mag={"cached-mag": str(faa)},
            direct_annotation=True,
            bakta_dir="",
            extension="faa",
        ),
        log=[str(log)],
    )

    runpy.run_path(str(cached_script), init_globals={"snakemake": fake_snakemake})

    staged = bins / "cached-mag.faa"
    assert staged.is_symlink()
    assert staged.resolve() == faa.resolve()
    assert str(cached_helper) not in sys.path


def test_two_bound_module_instances_are_deferred_and_isolated(tmp_path):
    configfile = tmp_path / "instances.yaml"
    configfile.write_text(
        yaml.safe_dump(
            {
                "first": _mock_config("first-instance"),
                "second": _mock_config("second-instance"),
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )
    harness = tmp_path / "two-instances.smk"
    harness.write_text(
        "PROVIDER_CALLS = []\n\n"
        "def first_sources(wildcards):\n"
        "    PROVIDER_CALLS.append('first')\n"
        "    return {'first-mag': 'upstream/first source/protein parent'}\n\n"
        "def first_taxonomy(wildcards):\n"
        "    return 'taxonomy/first.tsv'\n\n"
        "module first_modeling:\n"
        '    pathvars: results="results/first/{sample}"\n'
        f"    snakefile: {json.dumps(str(SNAKEFILE))}\n"
        '    config: config["first"]\n\n'
        "use rule * from first_modeling exclude all as first_*\n\n"
        "first_modeling.bind_protein_sources(\n"
        "    sources=first_sources, taxonomy=first_taxonomy\n"
        ")\n"
        "assert PROVIDER_CALLS == []\n\n"
        "module second_modeling:\n"
        '    pathvars: results="results/second/{sample}"\n'
        f"    snakefile: {json.dumps(str(SNAKEFILE))}\n"
        '    config: config["second"]\n\n'
        "use rule * from second_modeling exclude all as second_*\n\n"
        "second_modeling.bind_protein_sources(\n"
        "    sources={'second-mag': 'upstream/second source/protein parent'},\n"
        "    taxonomy='taxonomy/second.tsv',\n"
        ")\n"
        "assert PROVIDER_CALLS == []\n\n"
        "rule produce_protein_parent:\n"
        "    output:\n"
        "        directory('upstream/{instance} source/protein parent'),\n"
        "    params:\n"
        "        mag=lambda wc: wc.instance + '-mag',\n"
        "    shell:\n"
        "        \"mkdir -p {output:q}; printf '>{params.mag}\\nMPEPTIDE\\n' > {output:q}/annotation.faa\"\n\n"
        "rule produce_taxonomy:\n"
        "    output:\n"
        "        'taxonomy/{instance}.tsv',\n"
        "    params:\n"
        "        mag=lambda wc: wc.instance + '-mag',\n"
        "        domain=lambda wc: 'archaea' if wc.instance == 'first' else 'bacteria',\n"
        "    shell:\n"
        "        \"mkdir -p taxonomy; printf 'domain\\tuser_genome\\tclassification\\n{params.domain}\\t{params.mag}\\td__Mock\\n' > {output:q}\"\n\n"
        "rule all:\n"
        "    input:\n"
        "        first_modeling.FINAL_TARGETS,\n"
        "        second_modeling.FINAL_TARGETS,\n",
        encoding="utf-8",
    )

    listed = _run(tmp_path, harness, configfile, "--list-rules")
    assert listed.returncode == 0, listed.stdout + listed.stderr
    assert not (tmp_path / "upstream").exists()

    first = _run(tmp_path, harness, configfile, "--printshellcmds", "all")
    assert first.returncode == 0, first.stdout + first.stderr
    _assert_complete_mock_cohort(
        tmp_path, "first-instance", "first-mag", results_prefix="results/first"
    )
    _assert_complete_mock_cohort(
        tmp_path, "second-instance", "second-mag", results_prefix="results/second"
    )
    first_root = tmp_path / "results" / "first" / "first-instance"
    second_root = tmp_path / "results" / "second" / "second-instance"
    assert "first source/protein parent/annotation.faa" in (
        first_root / "qa" / "external_qa_summary.tsv"
    ).read_text(encoding="utf-8")
    assert "second source/protein parent/annotation.faa" not in (
        first_root / "qa" / "external_qa_summary.tsv"
    ).read_text(encoding="utf-8")
    assert "source=GTDB" in (
        first_root / "logs" / "carveme" / "first-mag.log"
    ).read_text(encoding="utf-8")
    assert "source=GTDB" in (
        second_root / "logs" / "carveme" / "second-mag.log"
    ).read_text(encoding="utf-8")
    rendered = first.stdout + first.stderr
    assert "first_sources" not in rendered

    cached = _run(tmp_path, harness, configfile, "all")
    assert cached.returncode == 0, cached.stdout + cached.stderr
    assert "Nothing to be done" in cached.stdout + cached.stderr


def test_bound_future_directories_support_automatic_checkm2(tmp_path):
    config = _mock_config("bound-auto-qa")
    config["qa"]["enabled"] = True
    config["qa"]["tsv"] = ""
    config["checkm2"]["database_path"] = str(tmp_path / "mock-checkm2-db")
    configfile = tmp_path / "auto-qa.yaml"
    configfile.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")

    fake_bin = _write_fake_checkm2(tmp_path)

    harness = tmp_path / "bound-auto-qa.smk"
    harness.write_text(
        "def sources(wildcards):\n"
        "    return {'future-mag': 'future producer/protein parent'}\n\n"
        "module modeling:\n"
        f"    snakefile: {json.dumps(str(SNAKEFILE))}\n"
        "    config: config\n\n"
        "use rule * from modeling exclude all as modeling_*\n\n"
        "modeling.bind_protein_sources(sources=sources)\n\n"
        "rule produce_future_parent:\n"
        "    output:\n"
        "        directory('future producer/protein parent'),\n"
        "    shell:\n"
        "        \"mkdir -p {output:q}; printf '>future-mag\\nMPEPTIDE\\n' > {output:q}/annotation.faa\"\n\n"
        "rule all:\n"
        "    input:\n"
        "        modeling.FINAL_TARGETS,\n",
        encoding="utf-8",
    )
    environment = {"PATH": f"{fake_bin}{os.pathsep}{os.environ['PATH']}"}

    result = _run(
        tmp_path,
        harness,
        configfile,
        "--printshellcmds",
        "all",
        environment=environment,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    _assert_complete_mock_cohort(tmp_path, "bound-auto-qa", "future-mag")
    staged = tmp_path / "results" / "bound-auto-qa" / "qa" / "checkm2" / "input"
    assert (staged / "future-mag.faa").is_symlink()
    assert (staged / "future-mag.faa").resolve() == (
        tmp_path / "future producer" / "protein parent" / "annotation.faa"
    ).resolve()
    rendered = result.stdout + result.stderr
    assert "localrule produce_future_parent:" in rendered
    assert "localrule modeling_stage_checkm2_faa:" in rendered
    assert "localrule modeling_checkm2:" in rendered
    assert "meew_metabolic_direct_inputs.py" in rendered
