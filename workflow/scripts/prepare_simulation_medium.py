"""Publish the one selected simulation medium and its source-model mapping."""

import json
import sys
import traceback
from pathlib import Path

sys.path.insert(0, str(snakemake.scriptdir))

from common.media import (
    aggregate_media,
    cobra_nutrient_mappings,
    model_ids_and_paths,
    read_medium_csv,
    write_mapping_tsv,
    write_medium_csv,
)


def _log(message):
    text = str(message)
    print(text, flush=True)
    log_path = Path(snakemake.log[0])
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a", encoding="utf-8") as handle:
        handle.write(text + "\n")


def main():
    source = str(snakemake.params.source)
    aggregation = str(snakemake.params.aggregation)
    tag = str(snakemake.params.tag)
    mock = bool(snakemake.params.mock)
    model_ids, model_paths = model_ids_and_paths(
        snakemake.input.models, snakemake.input.filtered_models
    )

    if source == "discover":
        recipes = list(snakemake.input.recipes)
        if len(recipes) != len(model_ids):
            raise ValueError(
                "Discovery requires exactly one valid recipe for every retained "
                f"model ({len(recipes)} recipe(s), {len(model_ids)} model(s))."
            )
        for model_id, recipe in zip(model_ids, recipes):
            if Path(recipe).parent.name != model_id:
                raise ValueError(
                    f"Recipe {recipe} does not belong to retained MAG {model_id!r}."
                )
        medium = aggregate_media(
            [read_medium_csv(recipe) for recipe in recipes], aggregation
        )
        selected_from = "individual minimum-total-import recipes"
    elif source == "predefined":
        predefined = list(snakemake.input.predefined)
        if len(predefined) != 1:
            raise ValueError(
                "Predefined mode requires exactly one declared medium CSV."
            )
        medium = read_medium_csv(predefined[0])
        selected_from = str(snakemake.params.predefined_label)
    else:
        raise ValueError(f"Unsupported simulation medium source: {source!r}")

    mapping_rows = []
    unavailable = {nutrient: [] for nutrient in medium}
    if not mock:
        from cobra.io import read_sbml_model

        raw_to_nutrient = {}
        for model_id, model_path in zip(model_ids, model_paths):
            model = read_sbml_model(model_path)
            mappings = cobra_nutrient_mappings(model, model_path, medium)
            for nutrient_id in medium:
                mapping = mappings.get(nutrient_id)
                if mapping is None:
                    unavailable[nutrient_id].append(model_id)
                    continue
                previous = raw_to_nutrient.setdefault(
                    mapping.raw_species_id, nutrient_id
                )
                if previous != nutrient_id:
                    raise ValueError(
                        f"Raw SBML species {mapping.raw_species_id!r} collides "
                        f"between canonical nutrients {previous!r} and "
                        f"{nutrient_id!r}."
                    )
                mapping_rows.append(
                    {
                        "mag": model_id,
                        "nutrient_id": nutrient_id,
                        "raw_species_id": mapping.raw_species_id,
                        "cobra_exchange_id": mapping.cobra_exchange_id,
                        "exchange_coefficient": mapping.exchange_coefficient,
                    }
                )
        unused = sorted(
            nutrient
            for nutrient, missing_models in unavailable.items()
            if len(missing_models) == len(model_ids)
        )
        if len(unused) == len(medium):
            raise ValueError(
                "None of the selected medium nutrients maps to an actual "
                "extracellular exchange in the retained model cohort."
            )
    else:
        unused = []
        unavailable = {nutrient: [] for nutrient in medium}

    write_medium_csv(snakemake.output.medium, medium)
    write_mapping_tsv(snakemake.output.mapping, mapping_rows)
    status = {
        "synthetic": mock,
        "tag": tag,
        "source": source,
        "aggregation": aggregation if source == "discover" else None,
        "predefined_medium": (
            str(snakemake.params.predefined_label) if source == "predefined" else None
        ),
        "selected_from": selected_from,
        "min_growth": float(snakemake.params.min_growth),
        "open_exchanges": bool(snakemake.params.open_exchanges),
        "exchange_bound": float(snakemake.params.exchange_bound),
        "model_ids": model_ids,
        "model_count": len(model_ids),
        "nutrient_count": len(medium),
        "unused_constituents": unused,
        "unavailable_models_by_nutrient": {
            nutrient: models for nutrient, models in unavailable.items() if models
        },
    }
    Path(snakemake.output.status).parent.mkdir(parents=True, exist_ok=True)
    Path(snakemake.output.status).write_text(
        json.dumps(status, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    _log(
        f"Published selected medium tag={tag!r}, source={source!r}, "
        f"nutrients={len(medium)}, unused={unused}."
    )


try:
    main()
except Exception:
    _log("Selected-medium preparation failed:\n" + traceback.format_exc())
    raise
