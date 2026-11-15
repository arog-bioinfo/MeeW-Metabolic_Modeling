"""Parse and normalize abundance tables for MICOM communities."""

import csv
import math
from dataclasses import dataclass
from pathlib import Path

COVERM_GENOME_COLUMN = "Genome"
COVERM_SAMPLE_COLUMN = "Sample"
COVERM_VALUE_COLUMN = "Relative Abundance (%)"
COVERM_DENSE_SUFFIX = f" {COVERM_VALUE_COLUMN}"
COVERM_UNMAPPED_GENOME = "unmapped"


@dataclass(frozen=True)
class NormalizedAbundances:
    """Abundances aligned to retained MICOM model order."""

    values: tuple[float, ...]
    extra_ids: tuple[str, ...]
    table_format: str


def _read_tsv(path):
    path = Path(path)
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t", strict=True)
        try:
            columns = next(reader)
        except StopIteration as error:
            raise ValueError(f"MICOM abundance TSV is empty: {path}") from error

        if not columns or any(column == "" for column in columns):
            raise ValueError("MICOM abundance TSV must have a non-empty header.")
        duplicate_columns = sorted(
            {column for column in columns if columns.count(column) > 1}
        )
        if duplicate_columns:
            raise ValueError(
                "MICOM abundance TSV contains duplicate column(s): "
                f"{duplicate_columns}"
            )

        rows = []
        for line_number, values in enumerate(reader, start=2):
            if len(values) != len(columns):
                raise ValueError(
                    f"MICOM abundance TSV line {line_number} has {len(values)} "
                    f"field(s); expected {len(columns)}."
                )
            rows.append(dict(zip(columns, values)))
    return columns, rows


def _select_columns(columns, rows):
    column_set = set(columns)
    if "abundance" in column_set and ({"id", "mag"} & column_set):
        return ("id" if "id" in column_set else "mag"), "abundance", "standard"

    sparse_required = {
        COVERM_SAMPLE_COLUMN,
        COVERM_GENOME_COLUMN,
        COVERM_VALUE_COLUMN,
    }
    if sparse_required.issubset(column_set):
        samples = {row[COVERM_SAMPLE_COLUMN] for row in rows}
        if "" in samples:
            raise ValueError(
                "Sparse CoverM abundance TSV contains an empty Sample value."
            )
        if len(samples) > 1:
            raise ValueError(
                "Sparse CoverM abundance TSV is ambiguous: expected exactly one "
                f"Sample value, found {sorted(samples)}."
            )
        return COVERM_GENOME_COLUMN, COVERM_VALUE_COLUMN, "coverm_sparse"

    if COVERM_GENOME_COLUMN in column_set:
        value_columns = [
            column
            for column in columns
            if column.endswith(COVERM_DENSE_SUFFIX)
            and len(column) > len(COVERM_DENSE_SUFFIX)
        ]
        if len(value_columns) != 1:
            raise ValueError(
                "Dense CoverM abundance TSV requires exactly one sample-prefixed "
                f"column ending in {COVERM_DENSE_SUFFIX!r}; found {value_columns}."
            )
        return COVERM_GENOME_COLUMN, value_columns[0], "coverm_dense"

    raise ValueError(
        "MICOM abundance TSV must contain either 'id' or 'mag' plus 'abundance', "
        "a dense CoverM table with 'Genome' and exactly one sample-prefixed "
        "' Relative Abundance (%)' column, or sparse CoverM columns 'Sample', "
        "'Genome', and 'Relative Abundance (%)'."
    )


def parse_and_normalize_abundances(path, retained_model_ids):
    """Return table abundances aligned and normalized to retained model IDs."""
    columns, rows = _read_tsv(path)
    id_column, value_column, table_format = _select_columns(columns, rows)

    selected_rows = [
        row
        for row in rows
        if not (
            table_format.startswith("coverm")
            and row[id_column] == COVERM_UNMAPPED_GENOME
        )
    ]
    empty_ids = [
        index for index, row in enumerate(selected_rows, start=2) if not row[id_column]
    ]
    if empty_ids:
        raise ValueError(
            f"MICOM abundance TSV contains an empty model ID on line(s): {empty_ids}"
        )

    ids = [row[id_column] for row in selected_rows]
    duplicate_ids = sorted({model_id for model_id in ids if ids.count(model_id) > 1})
    if duplicate_ids:
        raise ValueError(
            f"MICOM abundance TSV contains duplicate model ID(s): {duplicate_ids}"
        )

    abundance_by_id = {}
    invalid_ids = []
    negative_ids = []
    for row in selected_rows:
        model_id = row[id_column]
        try:
            value = float(row[value_column])
        except (TypeError, ValueError):
            invalid_ids.append(model_id)
            continue
        if not math.isfinite(value):
            invalid_ids.append(model_id)
        elif value < 0:
            negative_ids.append(model_id)
        else:
            abundance_by_id[model_id] = value

    if invalid_ids:
        raise ValueError(
            "MICOM abundance TSV contains NaN, infinite, or non-numeric abundance "
            f"for model ID(s): {sorted(set(invalid_ids))}"
        )
    if negative_ids:
        raise ValueError(
            "MICOM abundance TSV contains negative abundance for model ID(s): "
            f"{sorted(set(negative_ids))}"
        )

    retained_model_ids = list(retained_model_ids)
    missing_ids = sorted(set(retained_model_ids) - set(abundance_by_id))
    if missing_ids:
        raise ValueError(f"MICOM abundance TSV missing model ID(s): {missing_ids}")

    extra_ids = tuple(sorted(set(abundance_by_id) - set(retained_model_ids)))
    selected_values = [abundance_by_id[model_id] for model_id in retained_model_ids]
    try:
        total = math.fsum(selected_values)
    except OverflowError as error:
        raise ValueError(
            "MICOM abundance values for retained models must have a finite, "
            "positive total."
        ) from error
    if not math.isfinite(total) or total <= 0:
        raise ValueError(
            "MICOM abundance values for retained models must sum to a positive number."
        )

    return NormalizedAbundances(
        values=tuple(value / total for value in selected_values),
        extra_ids=extra_ids,
        table_format=table_format,
    )
