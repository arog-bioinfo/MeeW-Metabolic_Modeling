# Configuration

`config.yaml` controls the metabolic modeling workflow.
The normalized configuration is validated against
`workflow/schemas/config.schema.yaml` before the DAG is built. CheckM2 database
paths are required by the schema only when CheckM2 is enabled or selected as
the automatic QA provider; CPLEX configuration is conditional on
`solver: cplex`, and SCIP/mock runs do not require a CPLEX Studio path.

- `input.faa_samples`: user-supplied TSV with required `mag` and `path` columns plus an optional `universe` column. Trimmed nonblank values must be exactly `bacteria` or `archaea`; blank cells and two-column sheets default to `bacteria`. Invalid values report the MAG and line. See `faa_samples.example.tsv`.
- `input.faa_globs` / `input.faa_dir`: discover annotated `.faa` files. Nested globs such as `../MeeW-Annotation/results/bakta/*/*.faa` are supported.
- `input.bakta_dir` / `input.bakta_done`: discover Bakta batch `.faa` files recursively after an optional done marker. MAG IDs are `.faa` basenames/stems. Files matching `*hypothetical*.faa` are excluded case-insensitively. If no usable FAA files are discovered, the workflow fails clearly.
- `qa`: optional first-pass QA filter before CarveMe. Standalone Metabolic Modeling can consume an existing CheckM2-style report by configuring it as `qa.tsv` with `qa.enabled: true`. When imported by Main MeeW, Annotation owns CheckM2, and Main wires its report directly into this filter when both Annotation QA and Metabolic Modeling first-pass QA are enabled. For example, the standalone configuration is:

  ```yaml
  qa:
    enabled: true
    tsv: $PROJECT_DATA/checkm2/quality_report.tsv
  ```

  Default CheckM2-style columns are already `Name`, `Completeness`, and `Contamination`. Set `qa.enabled: true` with `qa.tsv: ""` to auto-run CheckM2 instead.

- `resources_root`: optional module resource root used by automatic CheckM2 path resolution. A nonblank value supplies `<resources_root>/checkm2_db`; blank or omitted preserves the previous ambient chain.
- `checkm2`: CheckM2 settings. CheckM2 auto-runs for first-pass QA when `qa.enabled: true` and `qa.tsv: ""`; setting `checkm2.enabled: true` also requests the CheckM2 report as a workflow target. Its database path uses explicit config, nonblank `resources_root`, `CHECKM2DB`, `MEEW_RESOURCES/checkm2_db`, `$XDG_DATA_HOME/meew/checkm2_db`, then `~/.local/share/meew/checkm2_db`. The database is not bundled and must be installed before CheckM2 runs.
- `solver` and `media`: solver used by model/community tools and CarveMe's gap-fill medium. CarveMe receives the selected `bacteria` or `archaea` through `--universe`, emits FBC2 explicitly, and SMETANA consumes models with FBC2 flavor explicitly.
- `cplex`: optional CPLEX Studio path and auto-upgrade toggle used only when `solver: cplex` to remove Community Edition limits inside Snakemake conda environments. An explicit `cplex.studio_dir` takes precedence over `CPLEX_STUDIO_DIR`; neither is required for another solver.
- `memote.min_total_score`: threshold used to select models for community analyses.
- `community.tools`: enable or disable SMETANA, MICOM, and SteadyCom outputs.
- `simulation_medium`: one medium shared by MICOM and SteadyCom. `source: discover` runs one independently rerunnable COBRA minimum-total-import LP per Memote-passing model using `min_growth`, `open_exchanges`, and `exchange_bound`, then applies unweighted `aggregation: sum` (default) or `max` across the complete retained cohort. `source: predefined` runs no individual optimization; `predefined_medium: M9` selects the bundled aerobic 20-component, oxygen-containing uptake-100 computational preset, while any other value is a user CSV path. Custom CSVs have exactly `nutrient_id,uptake`; nutrient IDs are COBRA extracellular metabolite IDs and uptake values are finite positive import-flux budgets, not concentrations. Discovery-only settings do not impose growth on predefined-mode individuals.
- `smetana.modes`: independent SMETANA mode targets. The default is `global: true` and `detailed: false`; when `community.tools.smetana` is enabled, at least one mode must be true. Disabling the tool targets neither output. Global and detailed retain `smetana/global.tsv` and `smetana/detailed.tsv`, with independent mode-named logs and benchmarks under `logs/smetana/`.
- `smetana.media_db` / `smetana.media`: optional custom-media selection. Supply both a TSV path and a nonempty list of names, or leave both empty to use SMETANA's default complete medium. The TSV must have a `medium` column containing every selected name. All selected names are passed comma-separated to one invocation with `--mediadb`; the workflow neither bundles a media database nor creates per-medium jobs.
- `micom.abundance_tsv`: optional declared material input. Existing TSVs may use `id` or `mag` plus `abundance`. Native CoverM 0.7 dense tables must have exact `Genome` plus exactly one sample-prefixed column ending in ` Relative Abundance (%)`; sparse tables must have exact `Sample`, `Genome`, and `Relative Abundance (%)` columns and exactly one sample value. The exact CoverM pseudo-genome `unmapped` is ignored. Duplicate IDs, invalid/negative values, missing retained model IDs, ambiguous CoverM columns, and a nonpositive selected total fail clearly. IDs removed before MICOM are reported and ignored, and only retained model abundances are normalized to sum to one. An empty setting preserves equal-abundance standalone behavior.
- `steadycom`: native ReFramed SteadyCom with `objective: community_growth`; variability uses objective fraction 0.9.

The workflow does not run old CheckM1 by default. Use an external QA TSV for first-pass filtering when available; otherwise CheckM2 is automatically used when first-pass QA is enabled and no TSV is provided.

When imported by Main MeeW, `bind_protein_sources` provides a private direct-directory adapter rather than any standalone input setting. The external-QA checkpoint receives matching deterministic MAG IDs plus the exact Annotation-owned Bakta directories, and every CarveMe job receives its corresponding directory as genome material. Both scripts validate only the fixed `annotation.faa` child after the parent directory dependency is available. If Annotation GTDB-Tk is enabled, Main additionally wires the imported rule's consolidated `classification.tsv` output into every selected CarveMe job. The module resolves an exact trimmed `user_genome`/`domain` row at job time; `bacteria` and `archaea` select the matching universe. Blank/unknown domains or no exact row produce a logged bacterial fallback, while a missing/unreadable promised file, malformed schema/row, or duplicate identifier fails. With GTDB-Tk disabled, there is no taxonomy dependency and the bacterial default is logged. This mode does not use `input.faa_samples`, globs, directory discovery, `input.bakta_dir`, or `input.bakta_done`; all standalone behaviors above remain unchanged when the direct adapter is absent.

After a full import and before DAG construction, call `metabolic_modeling.bind_protein_sources(sources=source_provider, taxonomy=taxonomy_provider)`. `sources` maps each MAG ID to its declared protein-producing parent directory; `taxonomy` is optional. Either may be a deferred callback, which remains module-controller state outside configuration and job parameters.

All configured paths and paths read from `input.faa_samples` expand `~` and
environment variables. Relative paths keep Snakemake working-directory
semantics, including `--directory .test` runs.

Real SMETANA jobs run one mode per job with verbose output. Each job copies its
selected models and optional media TSV into a unique directory below
Snakemake's automatic `tmpdir`, verifies only the expected mode file, atomically
publishes that file, and cleans the scratch directory.
