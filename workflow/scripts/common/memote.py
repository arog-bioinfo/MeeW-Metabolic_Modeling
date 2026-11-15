"""Shared parsing helpers for Memote score reports."""

import json
from pathlib import Path

MISSING_SCORE = "NA"


def memote_model_id(score_file):
    """Return the model identifier encoded in a Memote score path."""
    return Path(score_file).parent.name


def extract_section_scores(sections):
    """Convert Memote's section list into a mapping from section to score."""
    if not isinstance(sections, list):
        return {}
    return {
        section.get("section", ""): section.get("score", MISSING_SCORE)
        for section in sections
        if isinstance(section, dict)
    }


def average_annotation_score(section_scores):
    """Average numeric ``annotation_*`` subsection scores when present."""
    values = [
        value
        for section, value in section_scores.items()
        if section.startswith("annotation_") and isinstance(value, (int, float))
    ]
    return sum(values) / len(values) if values else MISSING_SCORE


def read_memote_report_row(score_file):
    """Parse one Memote score file into the shared aggregate-report columns."""
    with Path(score_file).open(encoding="utf-8") as handle:
        data = json.load(handle)

    score = data.get("score", {})
    section_scores = extract_section_scores(score.get("sections", []))
    row = {
        "mag": memote_model_id(score_file),
        "total_score": score.get("total_score", MISSING_SCORE),
        "consistency": section_scores.get("consistency", MISSING_SCORE),
        "annotation": average_annotation_score(section_scores),
    }
    return row, tuple(section_scores)
