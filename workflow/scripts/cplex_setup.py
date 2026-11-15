"""CPLEX license setup helpers for Snakemake script jobs."""

import contextlib
import fcntl
import hashlib
import os
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime


def _log(message, log_file=None):
    text = f"[cplex_setup {datetime.now().isoformat(timespec='seconds')}] {message}"
    print(text, file=sys.stderr, flush=True)
    if log_file:
        try:
            os.makedirs(os.path.dirname(log_file), exist_ok=True)
            with open(log_file, "a") as handle:
                handle.write(text + "\n")
        except OSError as error:
            print(
                f"[cplex_setup] Could not write to log file {log_file!r}: {error}",
                file=sys.stderr,
                flush=True,
            )


def _is_error_1016(error):
    return "1016" in str(error)


def _probe_cplex_size_limit():
    """Solve a tiny 1001-binary-variable model to detect CE size limits.

    Run the probe in a child process so an upgrade can replace CPLEX binaries
    before the workflow tool imports/loads CPLEX in this Python process.
    """
    probe = r"""
import sys
try:
    import cplex
except ImportError as error:
    print(error, file=sys.stderr)
    sys.exit(3)

problem = cplex.Cplex()
problem.set_log_stream(None)
problem.set_results_stream(None)
problem.set_warning_stream(None)
problem.set_error_stream(None)
problem.objective.set_sense(problem.objective.sense.minimize)
problem.variables.add(
    obj=[0.0] * 1001,
    lb=[0.0] * 1001,
    ub=[1.0] * 1001,
    types=[problem.variables.type.binary] * 1001,
    names=[f"x{i}" for i in range(1001)],
)
problem.solve()
"""
    result = subprocess.run(
        [sys.executable, "-c", probe], capture_output=True, text=True
    )
    if result.returncode == 0:
        return
    message = (result.stderr or result.stdout or "CPLEX probe failed").strip()
    if result.returncode == 3:
        raise RuntimeError(
            "solver is set to 'cplex', but the Python package 'cplex' is not "
            f"available in the active Snakemake conda environment: {message}"
        )
    raise RuntimeError(message)


def _upgrade_commands(studio_dir):
    yield [sys.executable, "-m", "docplex.cli.config", "--upgrade", studio_dir]
    docplex = shutil.which("docplex")
    if docplex:
        yield [docplex, "config", "--upgrade", studio_dir]


def _upgrade_lock_path():
    env_prefix = os.path.abspath(sys.prefix)
    lock_name = ".cplex_docplex_upgrade.lock"
    if os.access(env_prefix, os.W_OK):
        return os.path.join(env_prefix, lock_name)

    digest = hashlib.sha256(env_prefix.encode()).hexdigest()[:16]
    return os.path.join(tempfile.gettempdir(), f"cplex_docplex_upgrade_{digest}.lock")


@contextlib.contextmanager
def _cplex_upgrade_lock(log_file=None):
    lock_path = _upgrade_lock_path()
    _log(f"Waiting for CPLEX/docplex upgrade lock: {lock_path}", log_file)
    handle = None
    try:
        handle = open(lock_path, "a+")
        fcntl.flock(handle, fcntl.LOCK_EX)
    except OSError as error:
        if handle is not None:
            handle.close()
        raise RuntimeError(
            f"Unable to acquire CPLEX/docplex upgrade lock at {lock_path!r}: {error}"
        ) from error
    try:
        handle.seek(0)
        handle.truncate()
        handle.write(f"pid={os.getpid()}\nprefix={os.path.abspath(sys.prefix)}\n")
        handle.flush()
        _log(f"Acquired CPLEX/docplex upgrade lock: {lock_path}", log_file)
        yield
    finally:
        fcntl.flock(handle, fcntl.LOCK_UN)
        handle.close()
        _log(f"Released CPLEX/docplex upgrade lock: {lock_path}", log_file)


def _format_command_output(result):
    parts = []
    if result.stdout:
        parts.append(f"stdout:\n{result.stdout.strip()}")
    if result.stderr:
        parts.append(f"stderr:\n{result.stderr.strip()}")
    return "\n".join(parts) if parts else "no stdout/stderr captured"


def _run_docplex_upgrade(studio_dir, log_file=None):
    last_error = None
    for command in _upgrade_commands(studio_dir):
        _log(f"Running CPLEX/docplex upgrade command: {' '.join(command)}", log_file)
        try:
            result = subprocess.run(command, check=True, capture_output=True, text=True)
            if result.stdout or result.stderr:
                _log(
                    "CPLEX/docplex upgrade command output: "
                    f"{_format_command_output(result)}",
                    log_file,
                )
            return command
        except (OSError, subprocess.CalledProcessError) as error:
            last_error = error
            message = f"CPLEX/docplex upgrade command failed: {error}"
            if isinstance(error, subprocess.CalledProcessError):
                message += f"\n{_format_command_output(error)}"
            _log(message, log_file)
    raise RuntimeError(
        "Unable to run docplex config upgrade inside the active environment. "
        f"Tried python/docplex commands for CPLEX Studio path {studio_dir!r}. "
        f"Last error: {last_error}"
    )


def ensure_cplex_ready(solver, studio_dir, auto_upgrade=True, log_file=None):
    """Ensure CPLEX is usable beyond Community Edition limits when requested."""
    if solver != "cplex":
        return

    studio_dir = studio_dir or os.environ.get("CPLEX_STUDIO_DIR", "")
    if studio_dir:
        studio_dir = os.path.expanduser(os.path.expandvars(str(studio_dir)))

    try:
        _probe_cplex_size_limit()
        _log("CPLEX size-limit probe succeeded; no docplex upgrade needed.", log_file)
        return
    except Exception as error:
        if not _is_error_1016(error):
            raise RuntimeError(
                "CPLEX preflight failed before running the workflow tool. "
                f"Original error: {error}"
            ) from error
        first_error = error

    if not auto_upgrade:
        raise RuntimeError(
            "CPLEX appears to be using the Community Edition size limit "
            "(CPLEX Error 1016), and cplex.auto_upgrade is false. Run "
            f"'docplex config --upgrade {studio_dir or '<CPLEX_Studio_path>'}' "
            "inside the active Snakemake conda environment, or enable "
            "cplex.auto_upgrade."
        ) from first_error

    if not studio_dir:
        raise RuntimeError(
            "CPLEX appears to be using the Community Edition size limit "
            "(CPLEX Error 1016), but cplex.studio_dir is empty. Set "
            "cplex.studio_dir to your CPLEX Studio installation path or run "
            "docplex config --upgrade <CPLEX_Studio_path> manually inside the "
            "active Snakemake conda environment."
        ) from first_error

    with _cplex_upgrade_lock(log_file=log_file):
        try:
            _probe_cplex_size_limit()
            _log(
                "CPLEX size-limit probe succeeded after acquiring upgrade lock; "
                "docplex upgrade was already completed by another process.",
                log_file,
            )
            return
        except Exception as error:
            if not _is_error_1016(error):
                raise RuntimeError(
                    "CPLEX preflight failed after acquiring the docplex upgrade lock. "
                    f"Original error: {error}"
                ) from error

        command = _run_docplex_upgrade(studio_dir, log_file=log_file)

        try:
            _probe_cplex_size_limit()
            _log("CPLEX size-limit probe succeeded after docplex upgrade.", log_file)
        except Exception as error:
            raise RuntimeError(
                "CPLEX still reports a size-limit problem after running "
                f"'{' '.join(command)}'. Check that cplex.studio_dir points to a "
                "valid CPLEX Studio installation. Original post-upgrade error: "
                f"{error}"
            ) from error
