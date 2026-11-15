# --- Imports -----------------------------------------------------------------

import glob
import importlib.util
import math
import os
import re
from collections.abc import Mapping
from pathlib import Path

from snakemake.utils import validate


def _load_dag_helper(module_name, source):
    """Load one module-owned DAG helper without claiming a global package name."""
    spec = importlib.util.spec_from_file_location(module_name, source)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load Metabolic Modeling helper from {source}")
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    return helper


_metabolic_paths = _load_dag_helper(
    "_meew_metabolic_paths",
    Path(workflow.source_path("../scripts/meew_metabolic_paths.py")),
)
_smetana_helpers = _load_dag_helper(
    "_meew_metabolic_smetana",
    Path(workflow.source_path("../scripts/common/smetana.py")),
)
_path_collection_helpers = _load_dag_helper(
    "_meew_metabolic_path_collections",
    Path(workflow.source_path("../scripts/common/path_collections.py")),
)
# Materialize this job helper beside cached scripts for imported module execution.
_CARVEME_HELPER_SOURCE = Path(workflow.source_path("../scripts/common/carveme.py"))
_carveme_helpers = _load_dag_helper(
    "_meew_metabolic_carveme_inputs",
    _CARVEME_HELPER_SOURCE,
)

expand_path = _metabolic_paths.expand_path
resolve_cplex_path = _metabolic_paths.resolve_cplex_path
resolve_resource_path = _metabolic_paths.resolve_resource_path
normalize_smetana_config = _smetana_helpers.normalize_smetana_config
selected_smetana_modes = _smetana_helpers.selected_smetana_modes
validate_media_configuration = _smetana_helpers.validate_media_configuration
duplicate_path_stems = _path_collection_helpers.duplicate_path_stems
read_faa_samples = _carveme_helpers.read_faa_samples
_CARVEME_INDEX_HELPER_SOURCE = Path(
    workflow.source_path("../scripts/common/carveme_index.py")
)
_carveme_index_helpers = _load_dag_helper(
    "_meew_metabolic_carveme_index",
    _CARVEME_INDEX_HELPER_SOURCE,
)
_DIRECT_INPUTS_HELPER_SOURCE = Path(
    workflow.source_path("../scripts/meew_metabolic_direct_inputs.py")
)

input_config = config.setdefault("input", {})
for field in ["faa_samples", "faa_dir", "bakta_dir", "bakta_done"]:
    if field in input_config:
        input_config[field] = expand_path(input_config[field])
if "faa_globs" in input_config:
    configured_globs = input_config["faa_globs"]
    if configured_globs is None:
        configured_globs = []
    elif isinstance(configured_globs, str):
        configured_globs = [configured_globs]
    input_config["faa_globs"] = [expand_path(pattern) for pattern in configured_globs]

qa_path_config = config.get("qa", {})
if "tsv" in qa_path_config:
    qa_path_config["tsv"] = expand_path(qa_path_config["tsv"])

configured_resources_root = expand_path(config.get("resources_root", ""))
if "resources_root" in config:
    config["resources_root"] = configured_resources_root

checkm2_path_config = config.setdefault("checkm2", {})
checkm2_selected = bool(checkm2_path_config.get("enabled", False)) or (
    bool(qa_path_config.get("enabled", False)) and not qa_path_config.get("tsv", "")
)
if checkm2_selected:
    checkm2_path_config["database_path"] = resolve_resource_path(
        checkm2_path_config.get("database_path", ""),
        configured_resources_root,
        "CHECKM2DB",
        "checkm2_db",
    )
elif checkm2_path_config.get("database_path"):
    checkm2_path_config["database_path"] = expand_path(
        checkm2_path_config["database_path"]
    )

cplex_path_config = config.setdefault("cplex", {})
cplex_path_config["studio_dir"] = resolve_cplex_path(
    config.get("solver", "scip"), cplex_path_config.get("studio_dir", "")
)

if config.get("micom", {}).get("abundance_tsv"):
    config["micom"]["abundance_tsv"] = expand_path(config["micom"]["abundance_tsv"])

simulation_medium_config = config.setdefault("simulation_medium", {})
if not isinstance(simulation_medium_config, dict):
    raise ValueError("simulation_medium must be a configuration mapping.")
for key, value in {
    "source": "discover",
    "aggregation": "sum",
    "predefined_medium": "M9",
    "min_growth": 0.1,
    "open_exchanges": True,
    "exchange_bound": 1000,
}.items():
    simulation_medium_config.setdefault(key, value)
if (
    simulation_medium_config.get("source") == "predefined"
    and simulation_medium_config.get("predefined_medium") != "M9"
):
    simulation_medium_config["predefined_medium"] = expand_path(
        simulation_medium_config.get("predefined_medium", "")
    )

smetana_config = normalize_smetana_config(config)
smetana_config["media_db"] = expand_path(smetana_config.get("media_db", ""))

validate(config, schema="../schemas/config.schema.yaml")
for field in ("min_growth", "exchange_bound"):
    value = float(simulation_medium_config[field])
    if not math.isfinite(value):
        raise ValueError(f"simulation_medium.{field} must be finite, got {value!r}.")

SAMPLE = str(config.get("sample", "sample"))
RESULTS_ROOT = "results/{sample}"
LOGS_ROOT = "<results>/logs"

# --- Input helper functions --------------------------------------------------


_PROTEIN_SOURCES_PROVIDER = None
_TAXONOMY_PROVIDER = None
_BOUND_USER_MAG_PATTERN = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]{0,127}")
_BOUND_GENERATED_MAG_PATTERN = re.compile(
    r"[A-Za-z0-9][A-Za-z0-9._-]{0,127}_(?:prok|euk)_bin_[0-9]{3,}"
)
_BOUND_GENERATED_MAG_MAX_LENGTH = 192


def bind_protein_sources(*, sources, taxonomy=None):
    """Bind annotated-protein parent directories and optional taxonomy."""
    global _PROTEIN_SOURCES_PROVIDER, _TAXONOMY_PROVIDER
    if _PROTEIN_SOURCES_PROVIDER is not None:
        raise ValueError("Protein sources are already bound for this module instance.")
    if not callable(sources) and not isinstance(sources, Mapping):
        raise TypeError("Protein sources must be a mapping or deferred callback.")
    if (
        taxonomy is not None
        and not callable(taxonomy)
        and not isinstance(taxonomy, (str, os.PathLike))
    ):
        raise TypeError("Taxonomy must be a path, deferred callback, or None.")
    _PROTEIN_SOURCES_PROVIDER = sources
    _TAXONOMY_PROVIDER = taxonomy


def _protein_sources_are_bound():
    return _PROTEIN_SOURCES_PROVIDER is not None


def _provider_value(provider, wildcards):
    return provider(wildcards) if callable(provider) else provider


def _path_value(value, *, label):
    if not isinstance(value, (str, os.PathLike)):
        raise TypeError(f"{label} must be a string or path-like value.")
    path = os.fspath(value)
    if not isinstance(path, str):
        raise TypeError(f"{label} must resolve to a string path.")
    if not path:
        raise ValueError(f"{label} cannot be empty.")
    return path


def _validate_bound_mag_id(mag):
    if not isinstance(mag, str):
        raise TypeError("Bound protein source MAG IDs must be strings.")
    if _BOUND_USER_MAG_PATTERN.fullmatch(mag) is not None:
        return mag
    if (
        len(mag) <= _BOUND_GENERATED_MAG_MAX_LENGTH
        and _BOUND_GENERATED_MAG_PATTERN.fullmatch(mag) is not None
    ):
        return mag
    raise ValueError(f"Invalid bound protein source MAG ID: {mag!r}.")


def _bound_protein_sources(wildcards):
    supplied = _provider_value(_PROTEIN_SOURCES_PROVIDER, wildcards)
    if not isinstance(supplied, Mapping):
        raise ValueError("Bound protein source provider must return a mapping.")
    validated = {}
    for mag, directory in supplied.items():
        mag = _validate_bound_mag_id(mag)
        validated[mag] = _path_value(
            directory, label=f"Bound protein source directory for MAG {mag!r}"
        )
    return {mag: validated[mag] for mag in sorted(validated)}


def _bound_faa_by_mag(wildcards):
    return {
        mag: os.path.join(directory, "annotation.faa")
        for mag, directory in _bound_protein_sources(wildcards).items()
    }


def protein_source_directories_input(wildcards):
    if not _protein_sources_are_bound():
        return []
    return list(_bound_protein_sources(wildcards).values())


def direct_annotation_param(wildcards):
    return _protein_sources_are_bound()


def direct_inputs_helper_input(wildcards):
    if not _protein_sources_are_bound():
        return []
    return str(_DIRECT_INPUTS_HELPER_SOURCE)


def direct_annotation_mags_param(wildcards):
    if not _protein_sources_are_bound():
        return []
    return list(_bound_protein_sources(wildcards))


def standalone_bakta_dir_param(wildcards):
    if _protein_sources_are_bound():
        return ""
    return config.get("input", {}).get("bakta_dir", "")


def _as_list(value):
    if value is None:
        return []
    if isinstance(value, str):
        return [value]
    return list(value)


def _faa_globs():
    input_cfg = config.get("input", {})
    patterns = _as_list(input_cfg.get("faa_globs"))
    faa_dir = input_cfg.get("faa_dir")
    bakta_dir = input_cfg.get("bakta_dir")
    if faa_dir:
        patterns.append(os.path.join(faa_dir, "*.faa"))
    if bakta_dir:
        patterns.append(os.path.join(bakta_dir, "**", "*.faa"))
    return patterns or ["data/mags/*.faa"]


def _faa_samples_tsv():
    return config.get("input", {}).get("faa_samples", "")


def _bakta_done_marker():
    return config.get("input", {}).get("bakta_done", "")


def _is_non_hypothetical_faa(path):
    name = Path(path).name.lower()
    return name.endswith(".faa") and "hypothetical" not in name


def _read_faa_samples(tsv):
    return {
        mag: record["path"]
        for mag, record in read_faa_samples(tsv, expand_path=expand_path).items()
    }


# --- FAA discovery -----------------------------------------------------------


def _discover_faa():
    faa_samples = _faa_samples_tsv()
    if faa_samples:
        return _read_faa_samples(faa_samples)
    paths = []
    for pattern in _faa_globs():
        paths.extend(glob.glob(pattern, recursive=True))
    unique = sorted(
        {os.path.normpath(path) for path in paths if _is_non_hypothetical_faa(path)}
    )
    duplicates = duplicate_path_stems(unique)
    if duplicates:
        raise ValueError(f"Duplicate .faa basenames are not supported: {duplicates}")
    return {Path(path).stem: path for path in unique}


_faa_by_mag_cache = None


def _runtime_faa_by_mag():
    global _faa_by_mag_cache
    if _faa_by_mag_cache is None:
        _faa_by_mag_cache = _discover_faa()
        if not _faa_by_mag_cache:
            raise ValueError(
                "No .faa files found. Configure input.faa_samples, "
                "input.faa_globs, input.faa_dir, or input.bakta_dir in "
                "config/config.yaml."
            )
    return _faa_by_mag_cache


def _runtime_mags():
    return sorted(_runtime_faa_by_mag())


def faa_sample_sheet_input(wildcards):
    if _protein_sources_are_bound():
        return []
    return _faa_samples_tsv() or []


def gtdb_classification_input(wildcards):
    """Standalone discovery has no taxonomy material unless Main overrides it."""
    if _TAXONOMY_PROVIDER is None:
        return []
    supplied = _provider_value(_TAXONOMY_PROVIDER, wildcards)
    if supplied is None or supplied == "":
        return []
    if isinstance(supplied, (str, os.PathLike)):
        return _path_value(supplied, label="Bound taxonomy")
    try:
        values = list(supplied)
    except TypeError as error:
        raise TypeError(
            "Bound taxonomy provider must return one path or an empty value."
        ) from error
    if not values:
        return []
    if len(values) != 1:
        raise ValueError(
            "Bound taxonomy provider must return one path or an empty value."
        )
    return _path_value(values[0], label="Bound taxonomy")


def bakta_done_input(wildcards):
    if _protein_sources_are_bound():
        return []
    return _bakta_done_marker() or []


def checkm2_faa_inputs(wildcards):
    if _protein_sources_are_bound():
        return protein_source_directories_input(wildcards)
    if config.get("input", {}).get("bakta_dir"):
        return bakta_done_input(wildcards)
    if _faa_samples_tsv():
        return [_faa_samples_tsv()]
    return list(_runtime_faa_by_mag().values())


def checkm2_faa_by_mag_param(wildcards, input):
    if _protein_sources_are_bound():
        return _bound_faa_by_mag(wildcards)
    if config.get("input", {}).get("bakta_dir"):
        return {}
    return _runtime_faa_by_mag()


def external_qa_faa_by_mag_param(wildcards, input):
    if _protein_sources_are_bound():
        return {}
    if config.get("input", {}).get("bakta_dir"):
        return {}
    return _runtime_faa_by_mag()


def external_qa_mags_param(wildcards, input):
    if _protein_sources_are_bound():
        return direct_annotation_mags_param(wildcards)
    if config.get("input", {}).get("bakta_dir"):
        return []
    return _runtime_mags()


def faa_for_mag(wildcards):
    if _protein_sources_are_bound():
        sources = _bound_protein_sources(wildcards)
        if wildcards.mag not in sources:
            raise ValueError(
                f"MAG {wildcards.mag!r} is not present in bound protein sources."
            )
        return sources[wildcards.mag]
    return _runtime_faa_by_mag()[wildcards.mag]


def _format_rule_output(template, **wildcards):
    return str(template).format(**wildcards)


def _material_paths(value):
    if not value:
        return
    if isinstance(value, (str, Path)):
        yield Path(value)
        return
    for item in value:
        yield from _material_paths(item)


def _material_mtime_ns(path):
    """Mirror Snakemake's local mtime treatment for files and directories."""
    path = Path(path)
    details = path.lstat()
    mtimes = [details.st_mtime_ns]
    target = path.resolve() if path.is_symlink() else path
    target_details = target.stat()
    if path.is_symlink():
        mtimes.append(target_details.st_mtime_ns)
    if target.is_dir():
        timestamp = target / ".snakemake_timestamp"
        if timestamp.exists():
            mtimes.append(timestamp.stat().st_mtime_ns)
        else:
            mtimes.append(target_details.st_mtime_ns)
    return max(mtimes)


# --- QA configuration --------------------------------------------------------

qa_config = config.get("qa", {})
qa_enabled = bool(qa_config.get("enabled", False))
qa_tsv = qa_config.get("tsv", "")
qa_columns = qa_config.get("columns", {})
qa_id_col = qa_columns.get("id", "Name")
qa_completeness_col = qa_columns.get("completeness", "Completeness")
qa_contamination_col = qa_columns.get("contamination", "Contamination")
qa_min_completeness = float(qa_config.get("min_completeness", 50.0))
qa_max_contamination = float(qa_config.get("max_contamination", 10.0))
checkm2_config = config.get("checkm2", {})
checkm2_enabled = bool(checkm2_config.get("enabled", False))
qa_uses_checkm2 = qa_enabled and not qa_tsv
checkm2_report = "<results>/qa/checkm2/quality_report.tsv"


# --- External QA checkpoint helpers ------------------------------------------


def first_qa_filter_outputs(wildcards):
    return checkpoints.external_qa_filter.get(sample=SAMPLE).output


def first_qa_filtered_mags_file(wildcards):
    return first_qa_filter_outputs(wildcards).filtered


def first_qa_summary_file(wildcards):
    return first_qa_filter_outputs(wildcards).summary


def first_qa_passing_mags(wildcards):
    with open(first_qa_filtered_mags_file(wildcards), encoding="utf-8") as handle:
        mags = [line.strip() for line in handle if line.strip()]
    if qa_enabled and not mags:
        raise ValueError("First QA filtering removed all discovered .faa files.")
    return mags


def first_qa_passing_gem_files(wildcards):
    return [
        _format_rule_output(rules.carveme.output.model, sample=SAMPLE, mag=mag)
        for mag in first_qa_passing_mags(wildcards)
    ]


def first_qa_passing_memote_reports(wildcards):
    return [
        _format_rule_output(rules.memote.output.report, sample=SAMPLE, mag=mag)
        for mag in first_qa_passing_mags(wildcards)
    ]


def first_qa_passing_memote_scores(wildcards):
    return [
        _format_rule_output(rules.memote.output.score, sample=SAMPLE, mag=mag)
        for mag in first_qa_passing_mags(wildcards)
    ]


# --- Community configuration -------------------------------------------------

community_config = config.get("community", {})
community_enabled = bool(community_config.get("enabled", True))
community_tools = community_config.get("tools", {})
run_micom = community_enabled and bool(community_tools.get("micom", True))
run_steadycom = community_enabled and bool(community_tools.get("steadycom", True))
SMETANA_MODES = selected_smetana_modes(config)
validate_media_configuration(
    smetana_config.get("media_db", ""), smetana_config.get("media", [])
)
USE_CONDA_ENVS = not bool(config.get("execution", {}).get("mock", False))
CARVEME_INDEX_INVOCATION = _carveme_index_helpers.workflow_invocation_token(workflow)
CARVEME_INDEX_READINESS = (
    f".snakemake/meew-metabolic-modeling/carveme-index/"
    f"{CARVEME_INDEX_INVOCATION}.ready.json"
)
CARVEME_INDEX_LOG = (
    f".snakemake/meew-metabolic-modeling/carveme-index/"
    f"{CARVEME_INDEX_INVOCATION}.log"
)


def carveme_model_output(wildcards):
    return _format_rule_output(
        rules.carveme.output.model,
        sample=getattr(wildcards, "sample", SAMPLE),
        mag=wildcards.mag,
    )


def carveme_index_readiness_for_model(wildcards, *, model, materials):
    """Require preparation only for models needing mtime-based reconstruction."""
    if not USE_CONDA_ENVS:
        return []
    model_path = Path(model)
    try:
        model_details = model_path.stat()
    except OSError:
        return CARVEME_INDEX_READINESS
    if not model_path.is_file() or model_details.st_size == 0:
        return CARVEME_INDEX_READINESS
    for material in materials:
        for path in _material_paths(material):
            try:
                material_mtime_ns = _material_mtime_ns(path)
            except OSError:
                return CARVEME_INDEX_READINESS
            if material_mtime_ns > model_details.st_mtime_ns:
                return CARVEME_INDEX_READINESS
    return []


def carveme_index_readiness_input(wildcards):
    return carveme_index_readiness_for_model(
        wildcards,
        model=carveme_model_output(wildcards),
        materials=(
            faa_for_mag(wildcards),
            first_qa_filtered_mags_file(wildcards),
            faa_sample_sheet_input(wildcards),
            gtdb_classification_input(wildcards),
        ),
    )


# --- Memote checkpoint helper functions --------------------------------------


def memote_filter_outputs(wildcards):
    return checkpoints.filter_memote.get(sample=SAMPLE).output


def memote_passing_model_files(wildcards):
    filtered = memote_filter_outputs(wildcards).filtered
    with open(filtered, encoding="utf-8") as handle:
        mags = [line.strip() for line in handle if line.strip()]
    return [
        _format_rule_output(rules.carveme.output.model, sample=SAMPLE, mag=mag)
        for mag in mags
    ]


def memote_filtered_models_file(wildcards):
    return memote_filter_outputs(wildcards).filtered


def memote_quality_file(wildcards):
    return memote_filter_outputs(wildcards).quality


def smetana_media_db_input(wildcards):
    return smetana_config.get("media_db", "") or []


def micom_abundance_input(wildcards):
    return config.get("micom", {}).get("abundance_tsv", "") or []


# --- Shared simulation-medium helpers ----------------------------------------

_M9_MEDIUM_SOURCE = str(workflow.source_path("../resources/M9.csv"))
_workflow_configfiles = getattr(workflow, "configfiles", [])
SIMULATION_CONFIG_FILES = tuple(
    dict.fromkeys(str(path) for path in _workflow_configfiles)
)


def simulation_config_inputs(wildcards):
    """Declare user scientific config material for mtime-only profiles."""
    return list(SIMULATION_CONFIG_FILES)


simulation_medium_source = simulation_medium_config["source"]
simulation_medium_aggregation = simulation_medium_config["aggregation"]
simulation_predefined_value = simulation_medium_config["predefined_medium"]
simulation_predefined_label = (
    "M9" if simulation_predefined_value == "M9" else simulation_predefined_value
)
SIMULATION_MEDIUM_TAG = (
    ("summed_demand" if simulation_medium_aggregation == "sum" else "max_demand")
    if simulation_medium_source == "discover"
    else ("M9" if simulation_predefined_value == "M9" else "custom")
)


def selected_recipe_inputs(wildcards):
    if simulation_medium_source != "discover":
        return []
    return [
        _format_rule_output(
            rules.discover_individual_medium.output.medium,
            sample=SAMPLE,
            mag=Path(model).stem,
        )
        for model in memote_passing_model_files(wildcards)
    ]


def selected_predefined_medium_input(wildcards):
    if simulation_medium_source != "predefined":
        return []
    return [
        (
            _M9_MEDIUM_SOURCE
            if simulation_predefined_value == "M9"
            else simulation_predefined_value
        )
    ]


def selected_medium_output(name):
    return _format_rule_output(
        getattr(rules.prepare_simulation_medium.output, name), sample=SAMPLE
    )


def micom_output(name):
    return _format_rule_output(getattr(rules.micom.output, name), sample=SAMPLE)


def steadycom_output(name):
    return _format_rule_output(getattr(rules.steadycom.output, name), sample=SAMPLE)


# --- Workflow targets ---------------------------------------------------------


def final_targets(wildcards):
    # Resolve first QA before the downstream Memote checkpoint so a cached
    # upstream checkpoint is reevaluated before score inputs are materialized.
    first_qa_targets = [
        first_qa_summary_file(wildcards),
        first_qa_filtered_mags_file(wildcards),
        first_qa_passing_gem_files(wildcards),
        first_qa_passing_memote_reports(wildcards),
        first_qa_passing_memote_scores(wildcards),
    ]
    filter_outputs = memote_filter_outputs(wildcards)
    targets = [
        *first_qa_targets,
        _format_rule_output(rules.memote_summary.output.summary, sample=SAMPLE),
        filter_outputs.quality,
        filter_outputs.filtered,
        (
            [
                micom_output("community_growth"),
                micom_output("member_growth"),
                micom_output("fluxes"),
                micom_output("status"),
            ]
            if run_micom
            else []
        ),
        (
            _format_rule_output(
                rules.smetana_global.output.global_scores, sample=SAMPLE
            )
            if "global" in SMETANA_MODES
            else []
        ),
        (
            _format_rule_output(rules.smetana_detailed.output.detailed, sample=SAMPLE)
            if "detailed" in SMETANA_MODES
            else []
        ),
        (
            [
                steadycom_output("community_growth"),
                steadycom_output("abundances"),
                steadycom_output("exchanges"),
                steadycom_output("variability"),
                steadycom_output("status"),
            ]
            if run_steadycom
            else []
        ),
    ]
    if checkm2_enabled or qa_uses_checkm2:
        targets.append(_format_rule_output(rules.checkm2.output.report, sample=SAMPLE))
    return [item for target in targets for item in _as_list(target)]


FINAL_TARGETS = final_targets
