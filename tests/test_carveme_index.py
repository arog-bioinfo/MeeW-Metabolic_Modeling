import multiprocessing
import os
import stat
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

SCRIPTS = Path(__file__).resolve().parents[1] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

from common.carveme_index import (  # noqa: E402
    BundledIndex,
    CarveMeIndexError,
    ensure_bundled_index,
    verify_readiness,
    workflow_invocation_token,
    write_readiness,
)


def _identity(path):
    details = Path(path).stat()
    return {
        "path": str(Path(path).resolve()),
        "device": details.st_dev,
        "inode": details.st_ino,
        "size": details.st_size,
        "mtime_ns": details.st_mtime_ns,
    }


def _fixture_paths(root, *, database_root=None):
    root = Path(root)
    database_root = Path(database_root or root)
    return BundledIndex(
        root,
        root / "config.cfg",
        root / "bigg proteins.faa",
        database_root / "bigg proteins.dmnd",
    )


def _write_fake_diamond(tmp_path):
    executable = tmp_path / "fake bin" / "diamond"
    executable.parent.mkdir()
    executable.write_text(
        "#!/usr/bin/env python3\n"
        "import os, sys, time\n"
        "from pathlib import Path\n"
        "args = sys.argv[1:]\n"
        "record = Path(os.environ['FAKE_DIAMOND_RECORD'])\n"
        "if args[0] == 'makedb':\n"
        "    with record.open('a', encoding='utf-8') as handle:\n"
        "        handle.write('makedb\\n')\n"
        "    if Path(os.environ['FAKE_DIAMOND_FAIL']).exists():\n"
        "        raise SystemExit(29)\n"
        "    time.sleep(0.15)\n"
        "    output = Path(args[args.index('-d') + 1] + '.dmnd')\n"
        "    output.write_bytes(b'VALID-DIAMOND-INDEX\\n')\n"
        "elif args[0] == 'dbinfo':\n"
        "    database = Path(args[args.index('--db') + 1])\n"
        "    if not database.is_file() or not database.read_bytes().startswith(b'VALID'):\n"
        "        raise SystemExit(31)\n"
        "else:\n"
        "    raise SystemExit(37)\n",
        encoding="utf-8",
    )
    executable.chmod(executable.stat().st_mode | stat.S_IXUSR)
    return executable


def _ensure_worker(paths, diamond, log, queue):
    try:
        result = ensure_bundled_index(paths, diamond=diamond, log_path=log)
    except Exception as error:  # pragma: no cover - returned for parent assertion
        queue.put(("error", repr(error)))
    else:
        queue.put((result[0], result[3]))


def _write_source_material(paths):
    paths.config.write_text("native config\n", encoding="utf-8")
    paths.fasta.write_text(">protein\nMPEPTIDE\n", encoding="utf-8")


def test_workflow_invocation_token_reuses_only_one_workflow_object():
    tokens = iter(("first-workflow", "second-workflow"))
    first_workflow = SimpleNamespace()
    second_workflow = SimpleNamespace()

    first = workflow_invocation_token(first_workflow, token_factory=tokens.__next__)
    repeated = workflow_invocation_token(
        first_workflow,
        token_factory=lambda: pytest.fail("same workflow allocated a second token"),
    )
    second = workflow_invocation_token(second_workflow, token_factory=tokens.__next__)

    assert first == repeated == "first-workflow"
    assert second == "second-workflow"
    assert second != first


def test_multiprocess_cold_cache_builds_once_and_aliases_wait(tmp_path, monkeypatch):
    cache = tmp_path / "shared cache"
    cache.mkdir()
    alias = tmp_path / "cache alias"
    alias.symlink_to(cache, target_is_directory=True)
    paths = _fixture_paths(cache)
    alias_paths = _fixture_paths(cache, database_root=alias)
    _write_source_material(paths)
    diamond = _write_fake_diamond(tmp_path)
    record = tmp_path / "diamond calls.txt"
    monkeypatch.setenv("FAKE_DIAMOND_RECORD", str(record))
    monkeypatch.setenv("FAKE_DIAMOND_FAIL", str(tmp_path / "do-not-fail"))

    context = multiprocessing.get_context("fork")
    queue = context.Queue()
    processes = [
        context.Process(
            target=_ensure_worker,
            args=(
                paths if index % 2 == 0 else alias_paths,
                diamond,
                tmp_path / f"worker-{index}.log",
                queue,
            ),
        )
        for index in range(8)
    ]
    for process in processes:
        process.start()
    for process in processes:
        process.join(10)
        assert process.exitcode == 0

    results = [queue.get(timeout=2) for _ in processes]
    assert [status for status, _ in results].count("prepared") == 1
    assert [status for status, _ in results].count("reused") == 7
    assert record.read_text(encoding="utf-8").splitlines() == ["makedb"]
    assert paths.database.read_bytes() == b"VALID-DIAMOND-INDEX\n"
    assert not list(cache.glob(".meew-carveme-index-stage-*"))

    before = (paths.database.read_bytes(), paths.database.stat().st_mtime_ns)
    status, *_ = ensure_bundled_index(paths, diamond=diamond)
    assert status == "reused"
    assert (paths.database.read_bytes(), paths.database.stat().st_mtime_ns) == before
    assert record.read_text(encoding="utf-8").splitlines() == ["makedb"]


def test_failed_build_cleans_only_staging_and_lock_is_reusable(tmp_path, monkeypatch):
    paths = _fixture_paths(tmp_path)
    _write_source_material(paths)
    diamond = _write_fake_diamond(tmp_path)
    record = tmp_path / "diamond calls.txt"
    fail = tmp_path / "fail makedb"
    fail.write_text("fail\n", encoding="utf-8")
    monkeypatch.setenv("FAKE_DIAMOND_RECORD", str(record))
    monkeypatch.setenv("FAKE_DIAMOND_FAIL", str(fail))

    ready = tmp_path / "readiness.json"
    with pytest.raises(CarveMeIndexError, match="private staged index"):
        status, config, fasta, database = ensure_bundled_index(paths, diamond=diamond)
        write_readiness(
            ready,
            paths=paths,
            status=status,
            config_identity=config,
            fasta_identity=fasta,
            database_identity=database,
            initializer_log="",
        )
    assert not paths.database.exists()
    assert not ready.exists()
    assert not list(tmp_path.glob(".meew-carveme-index-stage-*"))

    fail.unlink()
    status, *_ = ensure_bundled_index(paths, diamond=diamond)
    assert status == "prepared"
    assert paths.database.read_bytes() == b"VALID-DIAMOND-INDEX\n"
    assert record.read_text(encoding="utf-8").splitlines() == [
        "makedb",
        "makedb",
    ]


def test_invalid_existing_index_is_preserved_and_never_rebuilt(tmp_path, monkeypatch):
    paths = _fixture_paths(tmp_path)
    _write_source_material(paths)
    paths.database.write_bytes(b"INCOMPLETE\n")
    os.utime(paths.database, ns=(1_700_000_000_000_000_000,) * 2)
    before = (paths.database.read_bytes(), paths.database.stat().st_mtime_ns)
    diamond = _write_fake_diamond(tmp_path)
    record = tmp_path / "diamond calls.txt"
    monkeypatch.setenv("FAKE_DIAMOND_RECORD", str(record))
    monkeypatch.setenv("FAKE_DIAMOND_FAIL", str(tmp_path / "do-not-fail"))

    with pytest.raises(
        CarveMeIndexError, match="pre-existing index was preserved unchanged"
    ):
        ensure_bundled_index(paths, diamond=diamond)

    assert (paths.database.read_bytes(), paths.database.stat().st_mtime_ns) == before
    assert not record.exists()
    assert not list(tmp_path.glob(".meew-carveme-index-stage-*"))


def test_readiness_rejects_missing_replaced_and_malformed_cache(tmp_path):
    paths = _fixture_paths(tmp_path)
    _write_source_material(paths)
    paths.database.write_bytes(b"VALID-DIAMOND-INDEX\n")
    ready = tmp_path / "ready.json"
    write_readiness(
        ready,
        paths=paths,
        status="reused",
        config_identity=_identity(paths.config),
        fasta_identity=_identity(paths.fasta),
        database_identity=_identity(paths.database),
        initializer_log="initializer.log",
    )
    assert verify_readiness(ready, paths=paths)["status"] == "reused"

    paths.database.unlink()
    with pytest.raises(CarveMeIndexError, match="missing or unreadable"):
        verify_readiness(ready, paths=paths)
    paths.database.write_bytes(b"VALID-REPLACEMENT\n")
    with pytest.raises(CarveMeIndexError, match="changed after initialization"):
        verify_readiness(ready, paths=paths)

    ready.write_text("", encoding="utf-8")
    with pytest.raises(CarveMeIndexError, match="missing or malformed"):
        verify_readiness(ready, paths=paths)
