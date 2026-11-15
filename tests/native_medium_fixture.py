"""Generate tiny solver-backed simulation-medium regression fixtures."""

import argparse
import json
import xml.etree.ElementTree as ET
from pathlib import Path

from cobra import Metabolite, Model, Reaction
from cobra.io import write_sbml_model


def _remove_oxygen_exchange(path, *, boundary):
    tree = ET.parse(path)
    root = tree.getroot()
    model_element = next(
        element for element in root.iter() if element.tag.rsplit("}", 1)[-1] == "model"
    )
    reaction_list = next(
        element
        for element in model_element
        if element.tag.rsplit("}", 1)[-1] == "listOfReactions"
    )
    oxygen_exchanges = [
        reaction
        for reaction in reaction_list
        if reaction.attrib.get("id") == "R_EX_o2_e"
    ]
    if oxygen_exchanges:
        assert len(oxygen_exchanges) == 1
        reaction_list.remove(oxygen_exchanges[0])
    oxygen = [
        species
        for species in root.iter()
        if species.tag.rsplit("}", 1)[-1] == "species"
        and species.attrib.get("id") == "M_o2_e"
    ]
    assert len(oxygen) == 1
    oxygen[0].set("boundaryCondition", "true" if boundary else "false")
    tree.write(path, encoding="utf-8", xml_declaration=True)


def _main_model(model_id, glucose_coefficient, output, *, block_secretion):
    model = Model(model_id)
    glc_e = Metabolite("glc__D_e", compartment="C_e")
    o2_e = Metabolite("o2_e", compartment="C_e")
    lac_e = Metabolite("lac__D_e", compartment="C_e")
    glc_c = Metabolite("glc__D_c", compartment="C_c")
    o2_c = Metabolite("o2_c", compartment="C_c")
    waste_c = Metabolite("waste_c", compartment="C_c")

    ex_glc = Reaction("EX_glc__D_e", lower_bound=-1000, upper_bound=1000)
    ex_glc.add_metabolites({glc_e: -1})
    ex_o2 = Reaction("EX_o2_e", lower_bound=-1000, upper_bound=1000)
    ex_o2.add_metabolites({o2_e: -1})
    ex_lac = Reaction(
        "EX_lac__D_e",
        lower_bound=-1000,
        upper_bound=0 if block_secretion else 1000,
    )
    ex_lac.add_metabolites({lac_e: -1})
    transport_glc = Reaction("T_GLC", lower_bound=0, upper_bound=1000)
    transport_glc.add_metabolites({glc_e: -1, glc_c: 1})
    transport_o2 = Reaction("T_O2", lower_bound=0, upper_bound=1000)
    transport_o2.add_metabolites({o2_e: -1, o2_c: 1})
    biomass = Reaction("BIOMASS", lower_bound=0, upper_bound=1000)
    biomass.add_metabolites({glc_c: -glucose_coefficient, o2_c: -1, waste_c: 1})
    waste = Reaction("WASTE", lower_bound=0, upper_bound=1000)
    waste.add_metabolites({waste_c: -1, lac_e: 1})
    model.add_reactions(
        [ex_glc, ex_o2, ex_lac, transport_glc, transport_o2, biomass, waste]
    )
    model.objective = biomass
    write_sbml_model(model, output)
    _remove_oxygen_exchange(output, boundary=True)


def _orphan_model(model_id, output):
    model = Model(model_id)
    glc_e = Metabolite("glc__D_e", compartment="C_e")
    o2_e = Metabolite("o2_e", compartment="C_e")
    ac_e = Metabolite("ac_e", compartment="C_e")
    lac_e = Metabolite("lac__D_e", compartment="C_e")
    glc_c = Metabolite("glc__D_c", compartment="C_c")
    ac_c = Metabolite("ac_c", compartment="C_c")
    precursor = Metabolite("precursor_c", compartment="C_c")
    waste_c = Metabolite("waste_c", compartment="C_c")

    ex_glc = Reaction("EX_glc__D_e", lower_bound=-1000, upper_bound=1000)
    ex_glc.add_metabolites({glc_e: -1})
    ex_ac = Reaction("EX_ac_e", lower_bound=-1000, upper_bound=1000)
    ex_ac.add_metabolites({ac_e: -1})
    ex_lac = Reaction("EX_lac__D_e", lower_bound=-1000, upper_bound=1000)
    ex_lac.add_metabolites({lac_e: -1})
    transport_glc = Reaction("T_GLC", lower_bound=0, upper_bound=1000)
    transport_glc.add_metabolites({glc_e: -1, glc_c: 1})
    transport_ac = Reaction("T_AC", lower_bound=0, upper_bound=1000)
    transport_ac.add_metabolites({ac_e: -1, ac_c: 1})
    use_glc = Reaction("USE_GLC", lower_bound=0, upper_bound=1000)
    use_glc.add_metabolites({glc_c: -1, precursor: 1})
    use_ac = Reaction("USE_AC", lower_bound=0, upper_bound=1000)
    use_ac.add_metabolites({ac_c: -1, precursor: 1})
    biomass = Reaction("BIOMASS", lower_bound=0, upper_bound=1000)
    biomass.add_metabolites({precursor: -1, waste_c: 1})
    waste = Reaction("WASTE", lower_bound=0, upper_bound=1000)
    waste.add_metabolites({waste_c: -1, lac_e: 1})
    model.add_reactions(
        [
            ex_glc,
            ex_ac,
            ex_lac,
            transport_glc,
            transport_ac,
            use_glc,
            use_ac,
            biomass,
            waste,
        ]
    )
    model.add_metabolites([o2_e])
    model.objective = biomass
    write_sbml_model(model, output)
    _remove_oxygen_exchange(output, boundary=False)


def _seed_contract(root, model_ids):
    data = root / "data"
    data.mkdir(parents=True, exist_ok=True)
    for model_id in model_ids:
        (data / f"{model_id}.faa").write_text(
            f">{model_id}\nMPEPTIDE\n", encoding="utf-8"
        )

    result = root / "results" / "native"
    qa = result / "qa"
    qa.mkdir(parents=True, exist_ok=True)
    (qa / "filtered_mags.txt").write_text(
        "".join(f"{model_id}\n" for model_id in model_ids), encoding="utf-8"
    )
    (qa / "external_qa_summary.tsv").write_text(
        "mag\tfaa\tcompleteness\tcontamination\tstatus\n"
        + "".join(
            f"{model_id}\tdata/{model_id}.faa\tNA\tNA\tPASS\n" for model_id in model_ids
        ),
        encoding="utf-8",
    )
    return result


def _seed_memote(result, model_ids):
    memote = result / "memote"
    for model_id in model_ids:
        folder = memote / model_id
        folder.mkdir(parents=True, exist_ok=True)
        (folder / "score.json").write_text(
            json.dumps({"score": {"total_score": 1.0}}) + "\n",
            encoding="utf-8",
        )
        (folder / "report.html").write_text("native fixture\n", encoding="utf-8")
    (memote / "memote_summary.tsv").write_text(
        "mag\ttotal_score\tconsistency\tannotation\n"
        + "".join(f"{model_id}\t1\t1\t1\n" for model_id in model_ids),
        encoding="utf-8",
    )
    (memote / "memote_quality.tsv").write_text(
        "mag\ttotal_score\tconsistency\tannotation\tstatus\n"
        + "".join(f"{model_id}\t1\t1\t1\tPASS\n" for model_id in model_ids),
        encoding="utf-8",
    )
    (memote / "filtered_models.txt").write_text(
        "".join(f"{model_id}\n" for model_id in model_ids), encoding="utf-8"
    )


def generate(root, scenario):
    if scenario == "main":
        model_ids = ("raw-A", "raw.B")
    elif scenario == "blocked":
        model_ids = ("blocked",)
    else:
        model_ids = ("orphan-one", "orphan-two")
    result = _seed_contract(root, model_ids)
    gems = result / "gems"
    gems.mkdir(parents=True, exist_ok=True)
    if scenario == "main":
        _main_model("raw-A", 1.0, gems / "raw-A.xml", block_secretion=False)
        _main_model("raw.B", 2.0, gems / "raw.B.xml", block_secretion=False)
        (root / "data" / "abundance.tsv").write_text(
            "id\tabundance\nraw-A\t1\nraw.B\t0\n", encoding="utf-8"
        )
    elif scenario == "blocked":
        _main_model("blocked", 1.0, gems / "blocked.xml", block_secretion=True)
    else:
        for model_id in model_ids:
            _orphan_model(model_id, gems / f"{model_id}.xml")
    _seed_memote(result, model_ids)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument(
        "--scenario", choices=("main", "blocked", "orphan"), required=True
    )
    args = parser.parse_args()
    generate(args.root, args.scenario)


if __name__ == "__main__":
    main()
