"""Small shared contracts for simulation media and source-model identities."""

import csv
import gzip
import math
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path

MEDIUM_COLUMNS = ("nutrient_id", "uptake")
MAPPING_COLUMNS = (
    "mag",
    "nutrient_id",
    "raw_species_id",
    "cobra_exchange_id",
    "exchange_coefficient",
)


@dataclass(frozen=True)
class RawSpecies:
    raw_id: str
    compartment: str
    boundary: bool


@dataclass(frozen=True)
class NutrientMapping:
    nutrient_id: str
    raw_species_id: str
    cobra_exchange_id: str
    exchange_coefficient: float


def read_medium_csv(path):
    """Read an exact two-column nutrient/import-budget CSV."""
    path = Path(path)
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, strict=True)
        if tuple(reader.fieldnames or ()) != MEDIUM_COLUMNS:
            raise ValueError(
                f"Medium CSV {path} must have exactly the columns "
                f"{list(MEDIUM_COLUMNS)} in that order."
            )
        medium = {}
        for line_number, row in enumerate(reader, start=2):
            if None in row:
                raise ValueError(
                    f"Medium CSV {path} must have exactly two values on line "
                    f"{line_number}."
                )
            nutrient_id = row["nutrient_id"].strip()
            if not nutrient_id:
                raise ValueError(
                    f"Medium CSV {path} has an empty nutrient_id on line "
                    f"{line_number}."
                )
            if nutrient_id in medium:
                raise ValueError(
                    f"Medium CSV {path} has duplicate nutrient_id "
                    f"{nutrient_id!r} on line {line_number}."
                )
            try:
                uptake = float(row["uptake"])
            except (TypeError, ValueError) as error:
                raise ValueError(
                    f"Medium CSV {path} has a non-numeric uptake for "
                    f"{nutrient_id!r} on line {line_number}."
                ) from error
            if not math.isfinite(uptake) or uptake <= 0:
                raise ValueError(
                    f"Medium CSV {path} uptake for {nutrient_id!r} must be "
                    "finite and positive."
                )
            medium[nutrient_id] = uptake
    if not medium:
        raise ValueError(f"Medium CSV {path} contains no nutrients.")
    return medium


def write_medium_csv(path, medium):
    """Write a deterministic canonical medium CSV."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=MEDIUM_COLUMNS)
        writer.writeheader()
        for nutrient_id in sorted(medium):
            uptake = float(medium[nutrient_id])
            if not math.isfinite(uptake) or uptake <= 0:
                raise ValueError(
                    f"Cannot write non-finite or non-positive uptake for "
                    f"{nutrient_id!r}."
                )
            writer.writerow({"nutrient_id": nutrient_id, "uptake": uptake})


def aggregate_media(media, method):
    """Aggregate complete per-model recipes by unweighted sum or maximum."""
    if method not in {"sum", "max"}:
        raise ValueError(f"Unsupported medium aggregation: {method!r}")
    media = list(media)
    if not media:
        raise ValueError("No individual medium recipes were supplied for aggregation.")
    nutrients = sorted({nutrient for recipe in media for nutrient in recipe})
    if method == "sum":
        return {
            nutrient: math.fsum(recipe.get(nutrient, 0.0) for recipe in media)
            for nutrient in nutrients
        }
    return {
        nutrient: max(recipe.get(nutrient, 0.0) for recipe in media)
        for nutrient in nutrients
    }


def model_ids_and_paths(models, filtered_models_file):
    """Return the complete Memote-retained cohort in checkpoint-file order."""
    model_paths = [models] if isinstance(models, (str, Path)) else list(models)
    by_id = {}
    for model_path in model_paths:
        model_id = Path(model_path).stem
        if model_id in by_id:
            raise ValueError(
                f"Duplicate Memote-passing model filename identity: {model_id!r}."
            )
        by_id[model_id] = str(model_path)

    retained = []
    seen = set()
    with Path(filtered_models_file).open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            model_id = line.strip()
            if not model_id:
                continue
            if model_id in seen:
                raise ValueError(
                    f"Memote filtered-model list contains duplicate MAG "
                    f"{model_id!r} on line {line_number}."
                )
            seen.add(model_id)
            retained.append(model_id)

    if not retained:
        raise ValueError("No Memote-passing models are available for simulation.")
    missing = sorted(set(retained) - set(by_id))
    unexpected = sorted(set(by_id) - set(retained))
    if missing or unexpected:
        raise ValueError(
            "The declared model cohort does not exactly match the Memote retained "
            f"cohort (missing model files={missing}, unexpected model files={unexpected})."
        )
    return retained, [by_id[model_id] for model_id in retained]


def backend_ids(model_ids):
    """Build an explicit reversible raw-MAG/backend identity mapping."""
    model_ids = list(model_ids)
    if len(model_ids) != len(set(model_ids)):
        raise ValueError("Raw MAG identifiers must be unique.")
    raw_to_backend = {
        model_id: f"member_{index:04d}"
        for index, model_id in enumerate(model_ids, start=1)
    }
    if len(set(raw_to_backend.values())) != len(raw_to_backend):
        raise ValueError("Generated community backend identifiers collide.")
    return raw_to_backend


def _open_xml(path):
    path = Path(path)
    return gzip.open(path, "rb") if path.suffix == ".gz" else path.open("rb")


def raw_sbml_species(path):
    """Read source SBML species IDs without inventing exchange identities."""
    with _open_xml(path) as handle:
        root = ET.parse(handle).getroot()
    species = []
    seen = set()
    for element in root.iter():
        if element.tag.rsplit("}", 1)[-1] != "species":
            continue
        raw_id = element.attrib.get("id", "")
        if not raw_id:
            raise ValueError(f"SBML model {path} contains a species without an ID.")
        if raw_id in seen:
            raise ValueError(f"SBML model {path} has duplicate species ID {raw_id!r}.")
        seen.add(raw_id)
        species.append(
            RawSpecies(
                raw_id=raw_id,
                compartment=element.attrib.get("compartment", ""),
                boundary=element.attrib.get("boundaryCondition", "false").lower()
                == "true",
            )
        )
    if not species:
        raise ValueError(f"SBML model {path} contains no species.")
    return species


def cobra_nutrient_mappings(model, model_path, nutrient_ids):
    """Map canonical nutrient IDs through actual COBRA exchanges and raw species."""
    from cobra.io import sbml as cobra_sbml
    from cobra.medium import find_external_compartment

    external = find_external_compartment(model)
    raw_by_cobra = {}
    for species in raw_sbml_species(model_path):
        cobra_id = cobra_sbml.F_REPLACE[cobra_sbml.F_SPECIE](species.raw_id)
        raw_by_cobra.setdefault(cobra_id, []).append(species)

    mappings = {}
    for nutrient_id in nutrient_ids:
        metabolites = [
            metabolite
            for metabolite in model.metabolites
            if metabolite.id == nutrient_id and metabolite.compartment == external
        ]
        if not metabolites:
            continue
        if len(metabolites) != 1:
            raise ValueError(
                f"Model {model_path} has ambiguous extracellular metabolite "
                f"identity for nutrient {nutrient_id!r}."
            )
        metabolite = metabolites[0]
        exchanges = [
            reaction
            for reaction in model.exchanges
            if metabolite in reaction.metabolites
        ]
        if not exchanges:
            # A declared extracellular species is not usable as a nutrient
            # unless the loaded model exposes an actual exchange for it.
            continue
        if len(exchanges) > 1:
            raise ValueError(
                f"Model {model_path} maps nutrient {nutrient_id!r} to "
                f"{len(exchanges)} COBRA exchange reactions; expected exactly one."
            )
        raw_candidates = [
            species
            for species in raw_by_cobra.get(nutrient_id, [])
            if species.compartment == metabolite.compartment
        ]
        if len(raw_candidates) != 1:
            raise ValueError(
                f"Model {model_path} maps nutrient {nutrient_id!r} to "
                f"{len(raw_candidates)} raw SBML species; expected exactly one."
            )
        exchange = exchanges[0]
        if len(exchange.metabolites) != 1:
            raise ValueError(
                f"Exchange {exchange.id!r} for nutrient {nutrient_id!r} in "
                f"{model_path} is not a one-metabolite exchange."
            )
        coefficient = float(exchange.metabolites[metabolite])
        if not math.isclose(abs(coefficient), 1.0, rel_tol=0.0, abs_tol=1e-12):
            raise ValueError(
                f"Exchange {exchange.id!r} for nutrient {nutrient_id!r} in "
                f"{model_path} has unsupported stoichiometric coefficient "
                f"{coefficient}; only unit scaling is supported."
            )
        mappings[nutrient_id] = NutrientMapping(
            nutrient_id=nutrient_id,
            raw_species_id=raw_candidates[0].raw_id,
            cobra_exchange_id=exchange.id,
            exchange_coefficient=coefficient,
        )
    return mappings


def write_mapping_tsv(path, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=MAPPING_COLUMNS, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row[column] for column in MAPPING_COLUMNS})


def read_mapping_tsv(path):
    path = Path(path)
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t", strict=True)
        if tuple(reader.fieldnames or ()) != MAPPING_COLUMNS:
            raise ValueError(
                f"Medium mapping TSV {path} must have exactly the columns "
                f"{list(MAPPING_COLUMNS)} in that order."
            )
        rows = []
        seen = set()
        for line_number, row in enumerate(reader, start=2):
            if None in row:
                raise ValueError(
                    f"Medium mapping TSV {path} has extra values on line "
                    f"{line_number}."
                )
            key = (row["mag"], row["nutrient_id"])
            if not all(row[column] for column in MAPPING_COLUMNS):
                raise ValueError(
                    f"Medium mapping TSV {path} has an empty value on line "
                    f"{line_number}."
                )
            if key in seen:
                raise ValueError(
                    f"Medium mapping TSV {path} has duplicate MAG/nutrient "
                    f"mapping {key!r}."
                )
            seen.add(key)
            try:
                coefficient = float(row["exchange_coefficient"])
            except ValueError as error:
                raise ValueError(
                    f"Medium mapping TSV {path} has invalid exchange coefficient "
                    f"on line {line_number}."
                ) from error
            if not math.isclose(abs(coefficient), 1.0, rel_tol=0.0, abs_tol=1e-12):
                raise ValueError(
                    f"Medium mapping TSV {path} has unsupported non-unit exchange "
                    f"coefficient on line {line_number}."
                )
            row["exchange_coefficient"] = coefficient
            rows.append(row)
    return rows
