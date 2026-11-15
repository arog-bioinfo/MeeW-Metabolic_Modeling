import math
import runpy
import sys
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

from common.abundance import parse_and_normalize_abundances  # noqa: E402
from micom_script import _publish_results, _run_micom_analyses  # noqa: E402


def _write_table(tmp_path, text):
    path = tmp_path / "abundance.tsv"
    path.write_text(text, encoding="utf-8")
    return path


@pytest.mark.parametrize(
    ("text", "table_format"),
    [
        ("id\tabundance\na\t2\nb\t6\n", "standard"),
        ("mag\tabundance\na\t2\nb\t6\n", "standard"),
        (
            "Genome\tsample Relative Abundance (%)\n"
            "a\t20\n"
            "unmapped\t70\n"
            "b\t60\n",
            "coverm_dense",
        ),
        (
            "Sample\tGenome\tRelative Abundance (%)\n"
            "sample\ta\t20\n"
            "sample\tunmapped\t70\n"
            "sample\tb\t60\n",
            "coverm_sparse",
        ),
    ],
)
def test_supported_abundance_tables_normalize_retained_models(
    tmp_path, text, table_format
):
    parsed = parse_and_normalize_abundances(_write_table(tmp_path, text), ["b", "a"])

    assert parsed.table_format == table_format
    assert parsed.extra_ids == ()
    assert parsed.values == pytest.approx((0.75, 0.25))
    assert math.fsum(parsed.values) == pytest.approx(1.0)


def test_extra_ids_are_reported_but_excluded_from_normalization(tmp_path):
    path = _write_table(tmp_path, "id\tabundance\na\t2\nb\t6\nremoved\t1000\n")

    parsed = parse_and_normalize_abundances(path, ["a", "b"])

    assert parsed.extra_ids == ("removed",)
    assert parsed.values == pytest.approx((0.25, 0.75))


@pytest.mark.parametrize(
    ("text", "message"),
    [
        ("name\tvalue\na\t1\n", "must contain either 'id' or 'mag'"),
        (
            "Genome\ta Relative Abundance (%)\tb Relative Abundance (%)\n" "a\t1\t2\n",
            "requires exactly one sample-prefixed column",
        ),
        (
            "Sample\tGenome\tRelative Abundance (%)\n" "one\ta\t1\n" "two\tb\t2\n",
            "expected exactly one Sample value",
        ),
        ("id\tabundance\na\t1\na\t2\n", "duplicate model ID"),
        ("id\tabundance\na\tnot-a-number\n", "non-numeric abundance"),
        ("id\tabundance\na\tNaN\n", "NaN, infinite"),
        ("id\tabundance\na\tinf\n", "NaN, infinite"),
        ("id\tabundance\na\t-1\n", "negative abundance"),
        ("id\tabundance\na\t1\n", "missing model ID"),
        ("id\tabundance\na\t0\nb\t0\n", "sum to a positive number"),
    ],
)
def test_invalid_abundance_tables_fail_clearly(tmp_path, text, message):
    path = _write_table(tmp_path, text)
    retained = (
        ["a", "b"] if "missing model ID" in message or "positive" in message else ["a"]
    )

    with pytest.raises(ValueError, match=message):
        parse_and_normalize_abundances(path, retained)


class _Named(dict):
    __getattr__ = dict.__getitem__


def _run_micom_mock(tmp_path, abundance=None):
    model_paths = [tmp_path / "a.xml", tmp_path / "b.xml"]
    for model in model_paths:
        model.touch()
    filtered = tmp_path / "filtered.txt"
    filtered.write_text("a\nb\n", encoding="utf-8")
    medium = tmp_path / "medium.csv"
    medium.write_text("nutrient_id,uptake\nmock_nutrient_e,1\n", encoding="utf-8")
    mapping = tmp_path / "mapping.tsv"
    mapping.write_text(
        "mag\tnutrient_id\traw_species_id\tcobra_exchange_id\t"
        "exchange_coefficient\n",
        encoding="utf-8",
    )
    medium_status = tmp_path / "medium-status.json"
    medium_status.write_text(
        '{"tag": "summed_demand", "source": "discover"}\n', encoding="utf-8"
    )
    community_growth = tmp_path / "community-growth.tsv"
    member_growth = tmp_path / "member-growth.tsv"
    fluxes = tmp_path / "fluxes.tsv"
    status = tmp_path / "status.json"
    inputs = _Named(
        models=[str(path) for path in model_paths],
        filtered_models=str(filtered),
        memote_quality=str(tmp_path / "quality.tsv"),
        abundance=[] if abundance is None else str(abundance),
        medium=str(medium),
        medium_mapping=str(mapping),
        medium_status=str(medium_status),
    )
    fake_snakemake = type(
        "FakeSnakemake",
        (),
        {
            "scriptdir": str(SCRIPTS),
            "input": inputs,
            "output": _Named(
                community_growth=str(community_growth),
                member_growth=str(member_growth),
                fluxes=str(fluxes),
                status=str(status),
            ),
            "params": _Named(solver="scip", mock=True, tag="summed_demand"),
            "config": {},
            "threads": 1,
            "log": [str(tmp_path / "micom.log")],
        },
    )()

    runpy.run_path(
        str(SCRIPTS / "micom_script.py"),
        run_name="__main__",
        init_globals={"snakemake": fake_snakemake},
    )
    return pd.read_csv(member_growth, sep="\t")


def test_micom_mock_preserves_no_table_equal_abundance_behavior(tmp_path, capsys):
    taxonomy = _run_micom_mock(tmp_path)

    assert taxonomy["stage"].tolist() == ["fba", "fba", "pfba", "pfba"]
    assert taxonomy["mag"].tolist() == ["a", "b", "a", "b"]
    assert taxonomy["abundance"].tolist() == pytest.approx([0.5] * 4)
    assert taxonomy["synthetic"].tolist() == [True] * 4
    assert "assigning equal abundances" in capsys.readouterr().out


def test_micom_mock_reads_coverm_from_declared_input_and_reports_extra_ids(
    tmp_path, capsys
):
    abundance = _write_table(
        tmp_path,
        "Genome\tsample Relative Abundance (%)\n"
        "a\t10\n"
        "b\t30\n"
        "removed\t60\n"
        "unmapped\t5\n",
    )

    taxonomy = _run_micom_mock(tmp_path, abundance)

    assert taxonomy["abundance"].tolist() == pytest.approx([0.25, 0.75, 0.25, 0.75])
    assert "Ignoring abundance IDs without Memote-passing models: ['removed']" in (
        capsys.readouterr().out
    )


def _solution(status="optimal", growth=2.0, objective=10.0, members=True):
    member_table = pd.DataFrame(
        {"abundance": [0.4, 0.6], "growth_rate": [1.0, 2.0]},
        index=["member_0001", "member_0002"],
    )
    if not members:
        member_table = member_table.iloc[:1]
    return SimpleNamespace(
        status=status,
        growth_rate=growth,
        objective_value=objective,
        members=member_table,
        fluxes=pd.DataFrame({"R_EX_glc": [-1.0, 1.0]}, index=["medium", "member_0001"]),
    )


class _InjectedNativeSolution:
    def __init__(self, owner):
        self.owner = owner

    def get_status_string(self):
        return self.owner.native_message


class _InjectedCommunity:
    def __init__(self, outcomes):
        self.outcomes = iter(outcomes)
        self.native_message = None
        self.solver = SimpleNamespace(
            status=None,
            problem=SimpleNamespace(solution=_InjectedNativeSolution(self)),
        )

    def optimize(self, **kwargs):
        outcome, native_message = next(self.outcomes)
        self.native_message = native_message
        if isinstance(outcome, BaseException):
            raise outcome
        self.solver.status = outcome.status
        return outcome


BACKEND_TO_RAW = {"member_0001": "mag-a", "member_0002": "mag-b"}


def test_micom_successful_pfba_keeps_complete_numerical_results():
    community = _InjectedCommunity(
        [(_solution(), "optimal"), (_solution(growth=1.9, objective=8.0), "optimal")]
    )

    community_rows, members, fluxes, overall, stages = _run_micom_analyses(
        community, BACKEND_TO_RAW, log=lambda message: None
    )

    assert overall == "complete"
    assert [(row["stage"], row["growth"]) for row in community_rows] == [
        ("fba", 2.0),
        ("pfba", 1.9),
    ]
    assert {row["stage"] for row in members} == {"fba", "pfba"}
    assert {row["stage"] for row in fluxes} == {"fba", "pfba"}
    assert stages["fba"]["result_status"] == "complete"
    assert stages["pfba"]["result_status"] == "complete"


def test_injected_nonoptimal_pfba_discards_finite_looking_pal6_values():
    community = _InjectedCommunity(
        [
            (_solution(), "optimal"),
            (
                _solution(
                    status="infeasible",
                    growth=2.099996966660152,
                    objective=36378.11954831087,
                ),
                "infeasible",
            ),
        ]
    )
    warnings = []

    community_rows, members, fluxes, overall, stages = _run_micom_analyses(
        community, BACKEND_TO_RAW, log=warnings.append
    )

    assert overall == "partial"
    assert community_rows[0]["stage"] == "fba" and community_rows[0]["growth"] == 2.0
    assert community_rows[1] == {
        "stage": "pfba",
        "growth": None,
        "objective_value": None,
        "stage_status": "unavailable",
        "synthetic": False,
    }
    assert {row["stage"] for row in members} == {"fba"}
    assert {row["stage"] for row in fluxes} == {"fba"}
    assert stages["pfba"]["native_status"] == "infeasible"
    assert stages["pfba"]["exception_type"] == "RuntimeError"
    assert "validated FBA outputs are retained" in warnings[0]
    assert "not exported" in warnings[0]


@pytest.mark.parametrize("native_message", ["numerical", "optimal"])
def test_injected_pfba_exception_does_not_inherit_stale_fba_diagnostic(
    native_message,
):
    community = _InjectedCommunity(
        [
            (_solution(), "optimal"),
            (RuntimeError("failure before a new solve"), native_message),
        ]
    )

    rows, members, fluxes, overall, stages = _run_micom_analyses(
        community, BACKEND_TO_RAW, log=lambda message: None
    )

    assert overall == "partial"
    assert rows[0]["stage"] == "fba" and rows[1]["growth"] is None
    assert {row["stage"] for row in members + fluxes} == {"fba"}
    assert stages["pfba"]["native_status"] is None
    assert stages["pfba"]["native_message"] is None
    assert stages["pfba"]["exception_type"] == "RuntimeError"
    assert stages["pfba"]["exception_message"] == "failure before a new solve"


@pytest.mark.parametrize(
    "outcome",
    [
        _solution(status="infeasible"),
        _solution(growth=math.nan),
        _solution(members=False),
        RuntimeError("primary optimizer failed"),
    ],
    ids=["nonoptimal", "nonfinite", "cohort-invalid", "optimizer-exception"],
)
def test_injected_primary_fba_failures_remain_fatal(outcome):
    community = _InjectedCommunity([(outcome, "primary failure")])

    with pytest.raises(RuntimeError):
        _run_micom_analyses(community, BACKEND_TO_RAW, log=lambda message: None)


def test_partial_micom_export_uses_na_and_preserves_fba_tables(tmp_path):
    community = _InjectedCommunity(
        [(_solution(), "optimal"), (_solution(status="infeasible"), "infeasible")]
    )
    rows, members, fluxes, overall, stages = _run_micom_analyses(
        community, BACKEND_TO_RAW, log=lambda message: None
    )
    output = SimpleNamespace(
        community_growth=tmp_path / "growth.tsv",
        member_growth=tmp_path / "members.tsv",
        fluxes=tmp_path / "fluxes.tsv",
        status=tmp_path / "status.json",
    )

    _publish_results(
        output,
        rows,
        members,
        fluxes,
        {"analysis_status": overall, "stage_status": stages},
    )

    growth = pd.read_csv(output.community_growth, sep="\t", keep_default_na=False)
    assert growth.loc[growth.stage == "fba", "growth"].item() == "2.0"
    failed = growth.loc[growth.stage == "pfba"].iloc[0]
    assert failed.growth == "NA" and failed.objective_value == "NA"
    assert failed.stage_status == "unavailable"
    assert set(pd.read_csv(output.member_growth, sep="\t").stage) == {"fba"}
    assert set(pd.read_csv(output.fluxes, sep="\t").stage) == {"fba"}
