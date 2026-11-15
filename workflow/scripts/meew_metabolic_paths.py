"""Portable path resolution for MeeW Metabolic Modeling."""

import os
from pathlib import Path


def expand_path(value):
    """Expand environment variables and ``~`` while preserving relative paths."""
    if value is None or str(value).strip() == "":
        return ""
    return os.path.expanduser(os.path.expandvars(str(value)))


def resources_root(configured=""):
    """Return a configured root or the existing ambient/XDG default chain."""
    if configured is not None and str(configured).strip():
        return Path(expand_path(configured))
    meew_resources = os.environ.get("MEEW_RESOURCES", "")
    if meew_resources.strip():
        return Path(expand_path(meew_resources))
    xdg_data_home = os.environ.get("XDG_DATA_HOME", "")
    if xdg_data_home.strip():
        return Path(expand_path(xdg_data_home)) / "meew"
    return Path(expand_path("~/.local/share/meew"))


def resolve_resource_path(
    explicit, configured_root, tool_environment, resource_subpath
):
    """Resolve explicit, configured root, tool environment, then ambient roots."""
    if explicit is not None and str(explicit).strip():
        return expand_path(explicit)
    if configured_root is not None and str(configured_root).strip():
        return str(resources_root(configured_root) / resource_subpath)
    tool_path = os.environ.get(tool_environment, "")
    if tool_path.strip():
        return expand_path(tool_path)
    return str(resources_root() / resource_subpath)


def resolve_optional_software_path(explicit, environment_name):
    """Resolve an optional software location without providing a default."""
    if explicit is not None and str(explicit).strip():
        return expand_path(explicit)
    environment_path = os.environ.get(environment_name, "")
    return expand_path(environment_path) if environment_path.strip() else ""


def resolve_cplex_path(solver, explicit):
    """Resolve CPLEX Studio only when the CPLEX solver is selected."""
    if solver != "cplex":
        return expand_path(explicit) if explicit else ""
    return resolve_optional_software_path(explicit, "CPLEX_STUDIO_DIR")
