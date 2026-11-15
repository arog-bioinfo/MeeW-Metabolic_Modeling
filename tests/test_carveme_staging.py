import importlib
import inspect
import json
import os
import runpy
import stat
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest
import yaml

SCRIPTS = Path(__file__).resolve().parents[1] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

from common.carveme import staged_faa  # noqa: E402


def material_snapshot(*paths):
    return {path: (path.read_bytes(), path.stat().st_mtime_ns) for path in paths}


def _file_identity(path):
    details = path.stat()
    return {
        "path": str(path.resolve()),
        "device": details.st_dev,
        "inode": details.st_ino,
        "size": details.st_size,
        "mtime_ns": details.st_mtime_ns,
    }


def _write_fake_carveme_distribution(root):
    site = root / "fake site"
    package = site / "carveme"
    generated = package / "data" / "generated"
    generated.mkdir(parents=True)
    config = package / "config.cfg"
    config.write_text(
        "[generated]\n"
        "fasta_file = data/generated/bigg_proteins.faa\n"
        "diamond_db = data/generated/bigg_proteins.dmnd\n",
        encoding="utf-8",
    )
    fasta = generated / "bigg_proteins.faa"
    fasta.write_text(">bigg\nMPEPTIDE\n", encoding="utf-8")
    database = generated / "bigg_proteins.dmnd"
    database.write_bytes(b"VALID-DIAMOND-INDEX\n")
    metadata = site / "carveme-1.6.1.dist-info"
    metadata.mkdir()
    (metadata / "METADATA").write_text(
        "Metadata-Version: 2.1\nName: carveme\nVersion: 1.6.1\n",
        encoding="utf-8",
    )
    return site, package, config, fasta, database


def _write_fake_readiness(root, monkeypatch):
    site, package, config, fasta, database = _write_fake_carveme_distribution(root)
    monkeypatch.syspath_prepend(str(site))
    ready = root / "carveme-index.ready.json"
    ready.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "status": "reused",
                "package": str(package.resolve()),
                "initializer_log": str(root / "initializer.log"),
                "config": _file_identity(config),
                "fasta": _file_identity(fasta),
                "database": _file_identity(database),
            }
        )
        + "\n",
        encoding="utf-8",
    )
    return ready


def test_staged_faa_preserves_source_and_sibling_and_uses_unique_scratch(tmp_path):
    source_dir = tmp_path / "annotation inputs"
    source_dir.mkdir()
    source = source_dir / "proteins with spaces.faa"
    source_bytes = b">protein\nMPEPTIDE\n"
    source.write_bytes(source_bytes)
    sibling = source.with_suffix(".tsv")
    sibling_bytes = b"Bakta annotation\n"
    sibling.write_bytes(sibling_bytes)
    other_source_file = source_dir / "keep.txt"
    other_source_file.write_bytes(b"unrelated Annotation output\n")
    for index, path in enumerate((source, sibling, other_source_file), start=1):
        os.utime(path, ns=(1_700_000_000_000_000_000 + index,) * 2)
    before = material_snapshot(source, sibling, other_source_file)

    scratch_root = tmp_path / "scratch root"
    other_job = scratch_root / "other-job"
    other_job.mkdir(parents=True)
    scratch_sentinel = other_job / "keep.txt"
    scratch_sentinel.write_bytes(b"other job\n")

    with staged_faa(source, scratch_root) as first_staged:
        first_workdir = first_staged.parent
        assert first_workdir.parent == scratch_root.resolve()
        assert first_staged.name == source.name
        assert first_staged.read_bytes() == source_bytes
        assert not first_staged.is_symlink()
        assert not os.path.samefile(first_staged, source)
        first_staged.with_suffix(".tsv").write_bytes(b"CarveMe DIAMOND output\n")

        with staged_faa(source, scratch_root) as second_staged:
            assert second_staged.parent.parent == scratch_root.resolve()
            assert second_staged.parent != first_workdir
            assert first_staged.exists()
            assert second_staged.exists()
            assert not os.path.samefile(second_staged, source)

    assert material_snapshot(source, sibling, other_source_file) == before
    assert scratch_sentinel.read_bytes() == b"other job\n"
    assert list(scratch_root.iterdir()) == [other_job]


def test_staged_faa_cleans_scratch_after_failure(tmp_path):
    source = tmp_path / "input.faa"
    source.write_bytes(b">protein\nM\n")
    scratch_root = tmp_path / "scratch"
    scratch_root.mkdir()
    scratch_sentinel = scratch_root / "keep.txt"
    scratch_sentinel.write_bytes(b"parent content\n")
    before = material_snapshot(source, scratch_sentinel)

    with pytest.raises(RuntimeError, match="CarveMe failed"):
        with staged_faa(source, scratch_root) as staged:
            staged.with_suffix(".tsv").write_bytes(b"partial DIAMOND output\n")
            raise RuntimeError("CarveMe failed")

    assert material_snapshot(source, scratch_sentinel) == before
    assert list(scratch_root.iterdir()) == [scratch_sentinel]


def test_staged_faa_rejects_invalid_source_and_scratch_without_source_changes(
    tmp_path,
):
    source = tmp_path / "input.faa"
    source.write_bytes(b">protein\nM\n")
    before = material_snapshot(source)
    missing_scratch = tmp_path / "unused scratch"

    with pytest.raises(FileNotFoundError, match="not a regular file"):
        with staged_faa(tmp_path / "missing.faa", missing_scratch):
            pass
    assert not missing_scratch.exists()

    scratch_file = tmp_path / "scratch-is-a-file"
    scratch_file.write_bytes(b"keep scratch file\n")
    scratch_before = material_snapshot(scratch_file)
    with pytest.raises(FileExistsError):
        with staged_faa(source, scratch_file):
            pass

    assert material_snapshot(source) == before
    assert material_snapshot(scratch_file) == scratch_before


@pytest.mark.parametrize("direct_annotation", [False, True])
def test_real_carveme_wrapper_invokes_carve_only_with_staged_faa(
    tmp_path, monkeypatch, direct_annotation
):
    workdir = tmp_path / "working directory"
    workdir.mkdir()
    monkeypatch.chdir(workdir)
    source_dir = tmp_path / "annotation source"
    source_dir.mkdir()
    source = source_dir / (
        "annotation.faa" if direct_annotation else "protein input.faa"
    )
    source_bytes = b">protein\nMPEPTIDE\n"
    source.write_bytes(source_bytes)
    sibling = source.with_suffix(".tsv")
    sibling_bytes = b"Bakta annotation\n"
    sibling.write_bytes(sibling_bytes)
    source_sentinel = source_dir / "other-output.txt"
    source_sentinel.write_bytes(b"keep me\n")
    faa_samples = ""
    gtdb_classification = ""
    if direct_annotation:
        classification = tmp_path / "annotation source" / "classification.tsv"
        classification.write_text(
            "domain\tuser_genome\tclassification\n" "archaea\tmag-a\td__Archaea\n",
            encoding="utf-8",
        )
        gtdb_classification = str(classification)
    else:
        sheet = tmp_path / "faa samples.tsv"
        sheet.write_text(
            f"mag\tpath\tuniverse\nmag-a\t{source}\tarchaea\n",
            encoding="utf-8",
        )
        faa_samples = str(sheet)
    protected = [source, sibling, source_sentinel]
    if direct_annotation:
        protected.append(classification)
    else:
        protected.append(sheet)
    for index, path in enumerate(protected, start=1):
        os.utime(path, ns=(1_710_000_000_000_000_000 + index,) * 2)
    source_before = material_snapshot(*protected)
    scratch_root = tmp_path / "scratch root"
    scratch_root.mkdir()
    scratch_sentinel = scratch_root / "other-job.txt"
    scratch_sentinel.write_bytes(b"keep other job\n")
    model = Path("published models/model.xml")
    log = Path("logs/carveme.log")
    index_ready = _write_fake_readiness(tmp_path / "index fixture", monkeypatch)

    cplex_setup = importlib.import_module("cplex_setup")
    cplex_calls = []
    monkeypatch.setattr(
        cplex_setup,
        "ensure_cplex_ready",
        lambda *args, **kwargs: cplex_calls.append((args, kwargs)),
    )

    shell_module = importlib.import_module("snakemake.shell")
    shell_calls = []
    staged_paths = []

    def carve_like_shell(command):
        staged = Path(inspect.currentframe().f_back.f_locals["staged_genome"])
        staged_paths.append(staged)
        assert staged != source
        assert staged.parent.parent == scratch_root.resolve()
        assert staged.name == source.name
        assert staged.read_bytes() == source_bytes
        assert not staged.is_symlink()
        assert not os.path.samefile(staged, source)
        staged.with_suffix(".tsv").write_bytes(b"CarveMe DIAMOND output\n")
        model.parent.mkdir(parents=True, exist_ok=True)
        model.write_text("model\n", encoding="utf-8")
        shell_calls.append(command)

    monkeypatch.setattr(shell_module, "shell", carve_like_shell)

    fake_snakemake = SimpleNamespace(
        scriptdir=SCRIPTS,
        input=SimpleNamespace(
            genome=str(source_dir if direct_annotation else source),
            faa_samples=faa_samples,
            gtdb_classification=gtdb_classification,
            index_ready=str(index_ready),
        ),
        output=SimpleNamespace(model=str(model)),
        params={
            "solver": "scip",
            "media": "M9",
            "mock": False,
            "direct_annotation": direct_annotation,
        },
        config={"cplex": {}},
        log=[str(log)],
        wildcards=SimpleNamespace(mag="mag-a"),
        log_fmt_shell=lambda **kwargs: f"> {log} 2>&1",
        resources=SimpleNamespace(tmpdir=str(scratch_root)),
    )

    monkeypatch.setitem(sys.modules, "common", SimpleNamespace())
    monkeypatch.delitem(sys.modules, "common.carveme", raising=False)
    runpy.run_path(
        str(SCRIPTS / "carveme.py"), init_globals={"snakemake": fake_snakemake}
    )

    assert len(staged_paths) == 1
    assert not staged_paths[0].exists()
    assert list(scratch_root.iterdir()) == [scratch_sentinel]
    assert material_snapshot(*protected) == source_before
    assert model.read_text(encoding="utf-8") == "model\n"
    assert cplex_calls == [(("scip", ""), {"auto_upgrade": True, "log_file": str(log)})]
    assert shell_calls == [
        "carve {staged_genome:q} --output {model:q} --fbc2 --universe "
        "{universe:q} --solver "
        "{solver:q} --gapfill {media:q} {log}"
    ]
    log_text = log.read_text(encoding="utf-8")
    assert "MAG=mag-a universe=archaea" in log_text
    assert ("source=GTDB" if direct_annotation else "source=explicit TSV") in log_text


@pytest.mark.parametrize("output_kind", ["absent", "empty"])
def test_native_exit_zero_without_nonempty_model_fails_with_log_reference(
    tmp_path, monkeypatch, output_kind
):
    source = tmp_path / "source proteins.faa"
    source.write_text(">protein\nMPEPTIDE\n", encoding="utf-8")
    source_before = material_snapshot(source)
    scratch = tmp_path / "scratch"
    model = tmp_path / "models" / "model.xml"
    log = tmp_path / "logs" / "carveme.log"
    index_ready = _write_fake_readiness(tmp_path / "index fixture", monkeypatch)

    cplex_setup = importlib.import_module("cplex_setup")
    monkeypatch.setattr(cplex_setup, "ensure_cplex_ready", lambda *args, **kwargs: None)
    shell_module = importlib.import_module("snakemake.shell")

    def exit_zero_without_model(command):
        if output_kind == "empty":
            model.parent.mkdir(parents=True)
            model.touch()

    monkeypatch.setattr(shell_module, "shell", exit_zero_without_model)
    fake_snakemake = SimpleNamespace(
        scriptdir=SCRIPTS,
        input=SimpleNamespace(
            genome=str(source),
            faa_samples="",
            gtdb_classification="",
            index_ready=str(index_ready),
        ),
        output=SimpleNamespace(model=str(model)),
        params={"solver": "scip", "media": "M9", "mock": False},
        config={"cplex": {}},
        log=[str(log)],
        wildcards=SimpleNamespace(mag="missing-model"),
        resources=SimpleNamespace(tmpdir=str(scratch)),
    )

    with pytest.raises(RuntimeError, match=r"non-empty model.*carveme\.log"):
        runpy.run_path(
            str(SCRIPTS / "carveme.py"), init_globals={"snakemake": fake_snakemake}
        )

    assert material_snapshot(source) == source_before
    assert not list(scratch.iterdir())
    assert not model.exists() or model.stat().st_size == 0


@pytest.mark.parametrize("direct_annotation", [False, True])
def test_mock_wrapper_never_stages_or_invokes_carveme(
    tmp_path, monkeypatch, direct_annotation
):
    source_dir = tmp_path / "mock annotation"
    source_dir.mkdir()
    source = source_dir / "annotation.faa"
    source.write_bytes(b">mock\nM\n")
    sibling = source.with_suffix(".tsv")
    sibling.write_bytes(b"original annotation\n")
    before = material_snapshot(source, sibling)
    scratch_root = tmp_path / "must-not-exist"
    model = tmp_path / "mock model.xml"
    log = tmp_path / "mock log.txt"

    cplex_setup = importlib.import_module("cplex_setup")
    monkeypatch.setattr(
        cplex_setup,
        "ensure_cplex_ready",
        lambda *args, **kwargs: pytest.fail("mock path ran CPLEX setup"),
    )
    shell_module = importlib.import_module("snakemake.shell")
    monkeypatch.setattr(
        shell_module,
        "shell",
        lambda *args, **kwargs: pytest.fail("mock path invoked CarveMe"),
    )
    fake_snakemake = SimpleNamespace(
        scriptdir=SCRIPTS,
        input=SimpleNamespace(
            genome=str(source_dir if direct_annotation else source),
            faa_samples="",
            gtdb_classification="",
            index_ready="",
        ),
        output=SimpleNamespace(model=str(model)),
        params={"mock": True, "direct_annotation": direct_annotation},
        config={"cplex": {}},
        log=[str(log)],
        wildcards=SimpleNamespace(mag="mock-mag"),
        log_fmt_shell=lambda **kwargs: "",
        resources=SimpleNamespace(tmpdir=str(scratch_root)),
    )

    with pytest.raises(SystemExit) as exit_info:
        runpy.run_path(
            str(SCRIPTS / "carveme.py"), init_globals={"snakemake": fake_snakemake}
        )

    assert exit_info.value.code is None
    assert not scratch_root.exists()
    assert material_snapshot(source, sibling) == before
    assert model.read_text(encoding="utf-8").startswith("<?xml")
    assert str(source) in log.read_text(encoding="utf-8")
    assert "MAG=mock-mag universe=bacteria source=default" in log.read_text(
        encoding="utf-8"
    )


def test_mock_wrapper_rejects_invalid_explicit_universe(tmp_path):
    source = tmp_path / "input.faa"
    source.write_text(">protein\nM\n", encoding="utf-8")
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        f"mag\tpath\tuniverse\ninvalid-mag\t{source}\tarchaeal\n",
        encoding="utf-8",
    )
    model = tmp_path / "model.xml"
    fake_snakemake = SimpleNamespace(
        scriptdir=SCRIPTS,
        input=SimpleNamespace(
            genome=str(source),
            faa_samples=str(sheet),
            gtdb_classification="",
            index_ready="",
        ),
        output=SimpleNamespace(model=str(model)),
        params={"mock": True},
        config={"cplex": {}},
        log=[str(tmp_path / "carveme.log")],
        wildcards=SimpleNamespace(mag="invalid-mag"),
        resources=SimpleNamespace(tmpdir=str(tmp_path / "scratch")),
    )

    with pytest.raises(ValueError, match=r"invalid-mag.*line 2.*bacteria.*archaea"):
        runpy.run_path(
            str(SCRIPTS / "carveme.py"), init_globals={"snakemake": fake_snakemake}
        )
    assert not model.exists()


def test_imported_rule_loads_cached_helper_and_fake_carve_cannot_touch_annotation(
    tmp_path,
):
    workdir = tmp_path / "imported work"
    workdir.mkdir()
    bakta = workdir / "Annotation results" / "bakta"
    bakta.mkdir(parents=True)
    source = bakta / "annotation.faa"
    source.write_bytes(b">protein header preserved\nMPEPTIDE\n")
    annotation_tsv = bakta / "annotation.tsv"
    annotation_tsv.write_bytes(b"Bakta annotation bytes\n")
    sentinel = bakta / "other output.txt"
    sentinel.write_bytes(b"other Annotation bytes\n")
    for index, path in enumerate((source, annotation_tsv, sentinel), start=1):
        os.utime(path, ns=(1_720_000_000_000_000_000 + index,) * 2)
    before = material_snapshot(source, annotation_tsv, sentinel)

    qa = workdir / "qa.txt"
    qa.write_text("mag\n", encoding="utf-8")
    scratch_root = workdir / "scratch parent"
    scratch_root.mkdir()
    scratch_sentinel = scratch_root / "other-job.txt"
    scratch_sentinel.write_bytes(b"other job stays\n")
    record = workdir / "fake carve record.json"
    model = workdir / "published models" / "model.xml"
    preserved_model = workdir / "published models" / "preserved.xml"
    log = workdir / "logs" / "carveme.log"
    classification = workdir / "Annotation results" / "gtdbtk" / "classification.tsv"
    classification.parent.mkdir(parents=True)
    classification.write_text(
        "domain\tuser_genome\tclassification\n" "archaea\tmag\td__Archaea\n",
        encoding="utf-8",
    )
    os.utime(classification, ns=(1_720_000_000_000_000_010,) * 2)
    classification_before = material_snapshot(classification)
    preserved_model.parent.mkdir(parents=True)
    preserved_model.write_text("preserved model\n", encoding="utf-8")
    preserved_before = material_snapshot(preserved_model)

    fake_bin = workdir / "fake bin"
    fake_bin.mkdir()
    fake_carve = fake_bin / "carve"
    fake_carve.write_text(
        "#!/usr/bin/env python3\n"
        "import json, os, sys\n"
        "from pathlib import Path\n"
        "args = sys.argv[1:]\n"
        "query = Path(args[0])\n"
        "source = Path(os.environ['CARVEME_TEST_SOURCE'])\n"
        "query.with_suffix('.tsv').write_bytes(b'fake DIAMOND sidecar\\n')\n"
        "output = Path(args[args.index('--output') + 1])\n"
        "output.parent.mkdir(parents=True, exist_ok=True)\n"
        "output.write_text('fake model\\n', encoding='utf-8')\n"
        "record = {\n"
        "    'args': args,\n"
        "    'query': str(query),\n"
        "    'query_bytes': query.read_text(encoding='utf-8'),\n"
        "    'query_is_symlink': query.is_symlink(),\n"
        "    'same_file': os.path.samefile(query, source),\n"
        "}\n"
        "Path(os.environ['CARVEME_TEST_RECORD']).write_text(json.dumps(record))\n",
        encoding="utf-8",
    )
    fake_carve.chmod(fake_carve.stat().st_mode | stat.S_IXUSR)
    (fake_bin / "python3").symlink_to(sys.executable)
    diamond_calls = workdir / "diamond calls.txt"
    fake_diamond = fake_bin / "diamond"
    fake_diamond.write_text(
        "#!/usr/bin/env python3\n"
        "import os, sys\n"
        "from pathlib import Path\n"
        "args = sys.argv[1:]\n"
        "with Path(os.environ['CARVEME_DIAMOND_CALLS']).open('a') as handle:\n"
        "    handle.write(args[0] + '\\n')\n"
        "if args[0] == 'makedb':\n"
        "    Path(args[args.index('-d') + 1] + '.dmnd').write_bytes(b'VALID-DIAMOND-INDEX\\n')\n"
        "elif args[0] == 'dbinfo':\n"
        "    database = Path(args[args.index('--db') + 1])\n"
        "    if not database.read_bytes().startswith(b'VALID'):\n"
        "        raise SystemExit(62)\n"
        "else:\n"
        "    raise SystemExit(61)\n",
        encoding="utf-8",
    )
    fake_diamond.chmod(fake_diamond.stat().st_mode | stat.S_IXUSR)
    fake_site, _, _, _, fake_database = _write_fake_carveme_distribution(
        workdir / "installed"
    )
    fake_database.unlink()
    blocked_calls = workdir / "blocked calls.txt"
    for name in (
        "bakta",
        "curl",
        "wget",
        "conda",
        "mamba",
        "micromamba",
        "pip",
        "pip3",
        "docplex",
    ):
        blocker = fake_bin / name
        blocker.write_text(
            "#!/bin/sh\n"
            f"printf '%s\\n' {name!r} >> \"$CARVEME_BLOCKED_CALLS\"\n"
            "exit 97\n",
            encoding="utf-8",
        )
        blocker.chmod(blocker.stat().st_mode | stat.S_IXUSR)

    config = yaml.safe_load(
        (SCRIPTS.parents[1] / "config" / "config.yaml").read_text(encoding="utf-8")
    )
    config["solver"] = "scip"
    config["execution"]["mock"] = False
    config["community"] = {
        "enabled": False,
        "tools": {"smetana": False, "micom": False, "steadycom": False},
    }
    snakefile = SCRIPTS.parent / "Snakefile"
    harness = workdir / "imported.smk"
    harness.write_text(
        "MODULE_CONFIG = "
        + repr(config)
        + "\n\nmodule metabolic_modeling:\n"
        + f"    snakefile: {json.dumps(str(snakefile))}\n"
        + "    config: MODULE_CONFIG\n"
        + "    skip_validation: True\n\n"
        + "use rule initialize_carveme_index from metabolic_modeling as imported_initialize_carveme_index\n\n"
        + "def imported_index_readiness(wildcards):\n"
        + "    return metabolic_modeling.carveme_index_readiness_for_model(\n"
        + f"        wildcards, model={json.dumps(str(model))},\n"
        + f"        materials=({json.dumps(str(bakta))}, {json.dumps(str(qa))}, "
        + f"{json.dumps(str(classification))}),\n"
        + "    )\n\n"
        + "def preserved_index_readiness(wildcards):\n"
        + "    return metabolic_modeling.carveme_index_readiness_for_model(\n"
        + f"        wildcards, model={json.dumps(str(preserved_model))},\n"
        + f"        materials=({json.dumps(str(bakta))}, {json.dumps(str(qa))}, "
        + f"{json.dumps(str(classification))}),\n"
        + "    )\n\n"
        + "use rule carveme from metabolic_modeling as imported_carveme with:\n"
        + "    input:\n"
        + f"        genome={json.dumps(str(bakta))},\n"
        + f"        qa={json.dumps(str(qa))},\n"
        + "        faa_samples=[],\n"
        + f"        gtdb_classification={json.dumps(str(classification))},\n"
        + "        index_ready=imported_index_readiness,\n"
        + "    output:\n"
        + f"        model={json.dumps(str(model))},\n"
        + "    log:\n"
        + f"        {json.dumps(str(log))}\n"
        + "    params:\n"
        + "        solver='scip',\n"
        + "        media='M9',\n"
        + "        mock=False,\n"
        + "        mag='mag',\n"
        + "        direct_annotation=True,\n"
        + "    resources:\n"
        + f"        tmpdir={json.dumps(str(scratch_root))}\n\n"
        + "use rule carveme from metabolic_modeling as preserved_carveme with:\n"
        + "    input:\n"
        + f"        genome={json.dumps(str(bakta))},\n"
        + f"        qa={json.dumps(str(qa))},\n"
        + "        faa_samples=[],\n"
        + f"        gtdb_classification={json.dumps(str(classification))},\n"
        + "        index_ready=preserved_index_readiness,\n"
        + "    output:\n"
        + f"        model={json.dumps(str(preserved_model))},\n"
        + "    log:\n"
        + f"        {json.dumps(str(workdir / 'logs' / 'preserved.log'))}\n"
        + "    params:\n"
        + "        solver='scip', media='M9', mock=False, mag='mag', direct_annotation=True,\n"
        + "    resources:\n"
        + f"        tmpdir={json.dumps(str(scratch_root))}\n",
        encoding="utf-8",
    )

    environment = os.environ.copy()
    environment.update(
        PATH=f"{fake_bin}{os.pathsep}/bin:/usr/bin",
        CARVEME_TEST_SOURCE=str(source),
        CARVEME_TEST_RECORD=str(record),
        CARVEME_DIAMOND_CALLS=str(diamond_calls),
        CARVEME_BLOCKED_CALLS=str(blocked_calls),
        CONDA_OFFLINE="true",
        PIP_NO_INDEX="1",
        PYTHONDONTWRITEBYTECODE="1",
        XDG_CACHE_HOME=str(workdir / "xdg cache"),
        PYTHONPATH=(str(fake_site) + os.pathsep + environment.get("PYTHONPATH", "")),
    )
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "snakemake",
            "--snakefile",
            str(harness),
            "--directory",
            str(workdir),
            "--cores",
            "1",
            "--nolock",
            str(model),
            str(preserved_model),
        ],
        check=False,
        capture_output=True,
        text=True,
        env=environment,
    )

    assert result.returncode == 0, result.stdout + result.stderr
    invocation = json.loads(record.read_text(encoding="utf-8"))
    staged_query = Path(invocation["query"])
    assert staged_query.name == source.name
    assert staged_query.parent.parent == scratch_root
    assert invocation["query_bytes"] == source.read_text(encoding="utf-8")
    assert invocation["query_is_symlink"] is False
    assert invocation["same_file"] is False
    assert invocation["args"] == [
        str(staged_query),
        "--output",
        str(model),
        "--fbc2",
        "--universe",
        "archaea",
        "--solver",
        "scip",
        "--gapfill",
        "M9",
    ]
    assert not staged_query.exists()
    assert list(scratch_root.iterdir()) == [scratch_sentinel]
    assert model.read_text(encoding="utf-8") == "fake model\n"
    assert material_snapshot(preserved_model) == preserved_before
    assert log.is_file()
    assert not (tmp_path / "imported").exists()
    assert material_snapshot(source, annotation_tsv, sentinel) == before
    assert material_snapshot(classification) == classification_before
    log_text = log.read_text(encoding="utf-8")
    assert "MAG=mag universe=archaea source=GTDB" in log_text
    assert f"source_path={classification}" in log_text
    assert diamond_calls.read_text(encoding="utf-8").splitlines() == [
        "makedb",
        "dbinfo",
    ]
    rendered = result.stdout + result.stderr
    assert rendered.count("localrule imported_initialize_carveme_index:") == 1
    assert rendered.index("rule imported_initialize_carveme_index:") < rendered.index(
        "rule imported_carveme:"
    )
    assert "localrule preserved_carveme:" not in rendered
    assert not blocked_calls.exists()


def test_failed_single_initializer_blocks_all_imported_native_models(tmp_path):
    workdir = tmp_path / "failed initialization"
    workdir.mkdir()
    sources = []
    for name in ("first", "second"):
        source = workdir / "inputs" / f"{name}.faa"
        source.parent.mkdir(exist_ok=True)
        source.write_text(f">{name}\nMPEPTIDE\n", encoding="utf-8")
        sources.append(source)
    source_before = material_snapshot(*sources)
    qa = workdir / "qa.txt"
    qa.write_text("first\nsecond\n", encoding="utf-8")

    fake_bin = workdir / "fake bin"
    fake_bin.mkdir()
    (fake_bin / "python3").symlink_to(sys.executable)
    diamond_record = workdir / "diamond calls.txt"
    diamond = fake_bin / "diamond"
    diamond.write_text(
        "#!/usr/bin/env python3\n"
        "import os, sys\n"
        "from pathlib import Path\n"
        "args = sys.argv[1:]\n"
        "if args[0] == 'makedb':\n"
        "    with Path(os.environ['CARVEME_DIAMOND_CALLS']).open('a') as handle:\n"
        "        handle.write('makedb\\n')\n"
        "    raise SystemExit(29)\n"
        "raise SystemExit(31)\n",
        encoding="utf-8",
    )
    diamond.chmod(diamond.stat().st_mode | stat.S_IXUSR)
    carve_record = workdir / "carve calls.txt"
    carve = fake_bin / "carve"
    carve.write_text(
        "#!/bin/sh\n" "printf 'carve\\n' >> \"$CARVEME_CARVE_CALLS\"\n" "exit 0\n",
        encoding="utf-8",
    )
    carve.chmod(carve.stat().st_mode | stat.S_IXUSR)
    fake_site, _, _, _, database = _write_fake_carveme_distribution(
        workdir / "installed"
    )
    database.unlink()

    config = yaml.safe_load(
        (SCRIPTS.parents[1] / "config" / "config.yaml").read_text(encoding="utf-8")
    )
    config["solver"] = "scip"
    config["execution"]["mock"] = False
    config["community"] = {
        "enabled": False,
        "tools": {"smetana": False, "micom": False, "steadycom": False},
    }
    snakefile = SCRIPTS.parent / "Snakefile"
    models = [workdir / "models" / f"{name}.xml" for name in ("first", "second")]
    harness = workdir / "failed-init.smk"
    harness.write_text(
        "MODULE_CONFIG = "
        + repr(config)
        + "\n\nmodule metabolic_modeling:\n"
        + f"    snakefile: {json.dumps(str(snakefile))}\n"
        + "    config: MODULE_CONFIG\n"
        + "    skip_validation: True\n\n"
        + "use rule initialize_carveme_index from metabolic_modeling as imported_initialize_carveme_index\n\n"
        + "def first_index_readiness(wildcards):\n"
        + "    return metabolic_modeling.carveme_index_readiness_for_model(\n"
        + f"        wildcards, model={json.dumps(str(models[0]))},\n"
        + f"        materials=({json.dumps(str(sources[0]))}, {json.dumps(str(qa))}),\n"
        + "    )\n\n"
        + "def second_index_readiness(wildcards):\n"
        + "    return metabolic_modeling.carveme_index_readiness_for_model(\n"
        + f"        wildcards, model={json.dumps(str(models[1]))},\n"
        + f"        materials=({json.dumps(str(sources[1]))}, {json.dumps(str(qa))}),\n"
        + "    )\n\n"
        + "use rule carveme from metabolic_modeling as first_model with:\n"
        + "    input:\n"
        + f"        genome={json.dumps(str(sources[0]))},\n"
        + f"        qa={json.dumps(str(qa))},\n"
        + "        faa_samples=[],\n"
        + "        gtdb_classification=[],\n"
        + "        index_ready=first_index_readiness,\n"
        + "    output:\n"
        + f"        model={json.dumps(str(models[0]))},\n"
        + "    log:\n"
        + f"        {json.dumps(str(workdir / 'first.log'))}\n"
        + "    params:\n"
        + "        solver='scip', media='M9', mock=False, mag='first', direct_annotation=False,\n\n"
        + "use rule carveme from metabolic_modeling as second_model with:\n"
        + "    input:\n"
        + f"        genome={json.dumps(str(sources[1]))},\n"
        + f"        qa={json.dumps(str(qa))},\n"
        + "        faa_samples=[],\n"
        + "        gtdb_classification=[],\n"
        + "        index_ready=second_index_readiness,\n"
        + "    output:\n"
        + f"        model={json.dumps(str(models[1]))},\n"
        + "    log:\n"
        + f"        {json.dumps(str(workdir / 'second.log'))}\n"
        + "    params:\n"
        + "        solver='scip', media='M9', mock=False, mag='second', direct_annotation=False,\n",
        encoding="utf-8",
    )
    environment = os.environ.copy()
    environment.update(
        PATH=f"{fake_bin}{os.pathsep}/bin:/usr/bin",
        PYTHONPATH=(str(fake_site) + os.pathsep + environment.get("PYTHONPATH", "")),
        PYTHONDONTWRITEBYTECODE="1",
        CARVEME_DIAMOND_CALLS=str(diamond_record),
        CARVEME_CARVE_CALLS=str(carve_record),
        XDG_CACHE_HOME=str(workdir / "xdg cache"),
    )
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "snakemake",
            "--snakefile",
            str(harness),
            "--directory",
            str(workdir),
            "--cores",
            "3",
            "--nolock",
            *(str(model) for model in models),
        ],
        check=False,
        capture_output=True,
        text=True,
        env=environment,
    )

    assert result.returncode != 0
    rendered = result.stdout + result.stderr
    assert rendered.count("localrule imported_initialize_carveme_index:") == 1
    assert "rule first_model:" not in rendered
    assert "rule second_model:" not in rendered
    assert diamond_record.read_text(encoding="utf-8").splitlines() == ["makedb"]
    assert not carve_record.exists()
    assert not any(model.exists() for model in models)
    assert not list(
        (workdir / ".snakemake" / "meew-metabolic-modeling" / "carveme-index").glob(
            "*.ready.json"
        )
    )
    assert material_snapshot(*sources) == source_before
