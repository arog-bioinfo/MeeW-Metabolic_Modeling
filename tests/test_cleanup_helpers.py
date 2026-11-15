import json
import sys
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parents[1] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

from common.community import (  # noqa: E402
    model_id_from_path,
    select_memote_passing_models,
)
from common.memote import read_memote_report_row  # noqa: E402
from common.path_collections import duplicate_path_stems  # noqa: E402


def test_memote_report_row_preserves_metrics_and_section_order(tmp_path):
    score_file = tmp_path / "model-a" / "score.json"
    score_file.parent.mkdir()
    score_file.write_text(
        json.dumps(
            {
                "score": {
                    "total_score": 0.8,
                    "sections": [
                        {"section": "consistency", "score": 0.9},
                        {"section": "annotation_metabolites", "score": 0.5},
                        {"section": "annotation_reactions", "score": 1.0},
                        {"section": "annotation_ignored", "score": "NA"},
                    ],
                }
            }
        ),
        encoding="utf-8",
    )

    row, section_names = read_memote_report_row(score_file)

    assert row == {
        "mag": "model-a",
        "total_score": 0.8,
        "consistency": 0.9,
        "annotation": 0.75,
    }
    assert section_names == (
        "consistency",
        "annotation_metabolites",
        "annotation_reactions",
        "annotation_ignored",
    )


def test_memote_report_row_keeps_missing_value_behavior(tmp_path):
    score_file = tmp_path / "model-b" / "score.json"
    score_file.parent.mkdir()
    score_file.write_text(
        json.dumps({"score": {"sections": "not-a-list"}}), encoding="utf-8"
    )

    row, section_names = read_memote_report_row(score_file)

    assert row == {
        "mag": "model-b",
        "total_score": "NA",
        "consistency": "NA",
        "annotation": "NA",
    }
    assert section_names == ()


def test_memote_report_row_propagates_malformed_json(tmp_path):
    score_file = tmp_path / "model-c" / "score.json"
    score_file.parent.mkdir()
    score_file.write_text("{", encoding="utf-8")

    with pytest.raises(json.JSONDecodeError):
        read_memote_report_row(score_file)


def test_community_model_selection_preserves_model_input_order(tmp_path):
    filtered_models = tmp_path / "filtered_models.txt"
    filtered_models.write_text("model-two\nmodel-one\nmodel-one\n", encoding="utf-8")
    models = [
        tmp_path / "excluded.xml",
        tmp_path / "model-one.xml",
        tmp_path / "model-two.xml",
    ]

    selected = select_memote_passing_models(models, filtered_models)

    assert selected == models[1:]
    assert [model_id_from_path(path) for path in selected] == [
        "model-one",
        "model-two",
    ]
    assert select_memote_passing_models(str(models[1]), filtered_models) == [
        str(models[1])
    ]


def test_duplicate_path_stems_are_unique_and_sorted():
    paths = [
        "z/beta.faa",
        "a/alpha.faa",
        "b/beta.faa",
        "c/alpha.faa",
        "d/beta.faa",
        "e/unique.faa",
    ]

    assert duplicate_path_stems(paths) == ["alpha", "beta"]
