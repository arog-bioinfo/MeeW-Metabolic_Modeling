# Snakemake pipeline: MeeW-Metabolic_Modeling

[![Snakemake](https://img.shields.io/badge/Snakemake-%E2%89%A59.12.0-039475?logo=snakemake&logoColor=white)](https://snakemake.readthedocs.io/)
[![run with conda](https://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/projects/conda/en/latest/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

<!-- TEMPORARY CATALOG OPT-OUT: later replace every standalone "pipeline" with the standard term formed by "work" + "flow", then remove this comment. -->

A Snakemake pipeline for reconstructing, assessing, and simulating prokaryotic genome-scale metabolic models.

> [!WARNING]
> **Research preview:** This pipeline is under active development and is not a deployment-ready release. Validate its results and resource requirements for your use case.

## Usage

The current input scope is prokaryotic annotated protein FASTA (`.faa`) files supplied by TSV, glob, directory, or Bakta output discovery. The pipeline can filter inputs with an external QA table or CheckM2, builds CarveMe models, assesses them with Memote, and runs enabled SMETANA, MICOM, and SteadyCom community analyses.

Install Snakemake `>=9.12.0` and Conda, then replace the input paths and review the solver, QA, and community settings in [`config/config.yaml`](config/config.yaml). Automatic CheckM2 QA requires its database. CPLEX is the default external solver and has separate eligibility, installation, and licensing requirements; native MICOM and SteadyCom execution requires a working CPLEX installation.

CarveMe's gap-filling `media` setting is separate from `simulation_medium`, the one shared medium used by MICOM and SteadyCom. That shared medium may be discovered with unweighted `sum` aggregation (default) or `max`, or predefined as M9 or a custom CSV. SMETANA medium handling is independent.

Results below `results/{sample}/` include QA decisions, CarveMe SBML models, Memote reports, and enabled community outputs. Valid primary MICOM or SteadyCom results can be retained when supplementary pFBA or variability analysis fails. Unavailable numerical results are recorded as `NA`, while `partial` and `unavailable` are explicit result or stage statuses. These retained partial results must not be interpreted as complete supplementary analyses or biological validation. The bundled mock fixture checks pipeline behavior only, not biological conclusions or full native solver execution.

See the [configuration guide](config/README.md) for input discovery, QA, solver setup, media semantics, and output contracts. Run the commands below from the repository root.

## Deployment options

Conda is the supported software deployment method. After configuring `config/config.yaml` for real inputs and an eligible solver, inspect the planned jobs and then execute them with:

```bash
snakemake --cores 2 --dry-run
snakemake \
  --cores 2 \
  --software-deployment-method conda
```

Adjust cores and resource overrides for the selected analyses and execution environment. This pipeline does not supply a container deployment.

## Execution profiles

The repository supplies [profile guidance](profiles/README.md), but no configured pipeline profile. Add or select a user-maintained profile for executor and site-specific resource settings.

## License

The pipeline is distributed under the [MIT License](LICENSE).

## References

No module-specific DOI is declared. Cite the repository version used and the selected modeling and assessment tools, including CheckM2, CarveMe, Memote, SMETANA, MICOM, SteadyCom/reframed, and Snakemake:

> Köster, J., Mölder, F., Jablonski, K. P., Letcher, B., Hall, M. B., Tomkins-Tinch, C. H., Sochat, V., Forster, J., Lee, S., Twardziok, S. O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., & Nahnsen, S. _Sustainable data analysis with Snakemake_. F1000Research, 10:33, **2021**. https://doi.org/10.12688/f1000research.29032.2.
