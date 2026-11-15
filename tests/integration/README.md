# Real-tool integration test

This harness executes the existing workflow against two unchanged, version-pinned public *E. coli* K-12 laboratory reference protein FASTAs. The external-QA values are synthetic, test-only metadata. This is reproducible software execution validation, not biological benchmarking or interpretation.

## Prerequisites and resources

- Run from the repository root with network access, `md5sum`, Conda, and the `snakemake-dev` environment.
- The complete scenario needs a working licensed IBM CPLEX Studio installation. Keep `cplex.studio_dir` empty in the tracked config and supply `CPLEX_STUDIO_DIR` only to the execution process.
- CheckM2 is disabled and is not needed. Snakemake creates isolated tool environments under the ignored `.test/.snakemake` tree.
- Allow several gigabytes for Conda environments and tool databases/caches. Downloads, environment solving, two CarveMe reconstructions, Memote, and community tools can take tens of minutes to hours depending on cache, network, CPUs, and solver setup. The commands use at most four cores.

## Prepare version-pinned fixtures

The preparation script downloads each gzip archive atomically, verifies its pinned MD5 before decompression, and writes the gzip payload unchanged to a stable `.faa` filename. Existing archives and FASTAs are verified on rerun. A mismatch stops preparation; remove only the reported corrupt ignored artifact before retrying.

```bash
mkdir -p .test/integration/run-logs
python tests/integration/prepare_fixtures.py \
  2>&1 | tee .test/integration/run-logs/prepare.log
md5sum -c tests/integration/fixture_archives.md5 \
  2>&1 | tee .test/integration/run-logs/md5.log
```

The IDs are the unambiguous FASTA stems `ecoli_k12_mg1655_GCF_000005845.2` and `ecoli_k12_w3110_GCF_000010245.2`. URLs, accessions, and expected hashes are tracked in `fixture_manifest.tsv`.

## Targeted external-QA filter

`synthetic_one_pass_one_fail.tsv` is deliberately configured for MG1655 to pass and W3110 to fail. Request only the checkpoint's QA outputs; do not use the default or a modeling target in this scenario.

```bash
conda run -n snakemake-dev snakemake \
  --directory .test \
  --configfile .test/integration/config/filter_only.yaml \
  --cores 1 --sdm conda --printshellcmds --show-failed-logs \
  results/integration_filter_only/qa/external_qa_summary.tsv \
  results/integration_filter_only/qa/filtered_mags.txt \
  2>&1 | tee .test/integration/run-logs/filter.log

PYTHONDONTWRITEBYTECODE=1 python tests/integration/assert_integration.py filter \
  --results .test/results/integration_filter_only \
  2>&1 | tee .test/integration/run-logs/filter-assertions.log
```

## Complete all-pass real-tool scenario

`synthetic_all_pass.tsv` passes both fixtures. The config uses `execution.mock: false`, CPLEX, no CheckM2, a zero Memote selection floor so any valid numeric real score can proceed, and all existing community tools. Replace the example CPLEX path with the local installation:

```bash
CPLEX_STUDIO_DIR=/path/to/CPLEX_Studio2212 \
conda run -n snakemake-dev snakemake \
  --directory .test \
  --configfile .test/integration/config/all_pass.yaml \
  --cores 4 --sdm conda --printshellcmds --show-failed-logs \
  --rerun-incomplete \
  2>&1 | tee .test/integration/run-logs/full.log
```

Stop at and report the first existing workflow, environment, solver, or license failure; do not change scientific algorithms to make the run pass. After a failure, validate only artifacts from stages that completed:

```bash
PYTHONDONTWRITEBYTECODE=1 python tests/integration/assert_integration.py full \
  --results .test/results/integration_all_pass \
  2>&1 | tee .test/integration/run-logs/full-partial-assertions.log
```

After successful completion, require every contract:

```bash
PYTHONDONTWRITEBYTECODE=1 python tests/integration/assert_integration.py full \
  --results .test/results/integration_all_pass --require-complete \
  2>&1 | tee .test/integration/run-logs/full-assertions.log
```

The assertions reject mock markers and placeholder SBML; check exact QA IDs/statuses, parseable models, numeric Memote totals, and the default `summed_demand` native outputs. MICOM must provide validated FBA plus explicit pFBA data or unavailability diagnostics. SteadyCom must provide positive community growth, complete abundances, actual finite nonzero exchanges, and explicit 0.9 SteadyComVA/status reporting. Without `--require-complete`, the harness accepts correctly marked partial pFBA or variability data with literal `NA` values and diagnostics; `--require-complete` rejects either partial result instead of calling it complete.

## Cleanup

All downloaded or generated material is ignored. Remove it without touching tracked harness files:

```bash
rm -rf tests/integration/fixtures \
  .test/integration/run-logs \
  .test/results/integration_filter_only \
  .test/results/integration_all_pass \
  .test/.snakemake
```
