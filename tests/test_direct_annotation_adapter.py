import sys
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parents[1] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

from meew_metabolic_direct_inputs import (  # noqa: E402
    annotation_faa,
    annotation_faa_by_mag,
)


def test_direct_annotation_adapter_uses_only_fixed_child(tmp_path):
    bakta = tmp_path / "mag-a" / "bakta"
    nested = bakta / "nested"
    nested.mkdir(parents=True)
    (nested / "discovered.faa").write_text(">wrong\nM\n", encoding="utf-8")

    with pytest.raises(FileNotFoundError, match="fixed 'annotation.faa' child"):
        annotation_faa(bakta)

    expected = bakta / "annotation.faa"
    expected.write_text(">right\nM\n", encoding="utf-8")
    assert annotation_faa(bakta) == str(expected)


def test_direct_annotation_adapter_pairs_caller_ids_and_directories_in_order(tmp_path):
    directories = []
    for name in ("first-parent", "second-parent"):
        directory = tmp_path / name
        directory.mkdir()
        (directory / "annotation.faa").write_text(f">{name}\nM\n", encoding="utf-8")
        directories.append(directory)

    assert annotation_faa_by_mag(["mag-z", "mag-a"], directories) == {
        "mag-z": str(directories[0] / "annotation.faa"),
        "mag-a": str(directories[1] / "annotation.faa"),
    }


@pytest.mark.parametrize(
    ("mags", "directories", "message"),
    [
        (["one"], [], "1 MAG IDs and 0 Bakta directories"),
        (["one", "one"], ["a", "b"], "duplicate MAG IDs"),
    ],
)
def test_direct_annotation_adapter_rejects_mismatched_identity_contracts(
    mags, directories, message
):
    with pytest.raises(ValueError, match=message):
        annotation_faa_by_mag(mags, directories)
