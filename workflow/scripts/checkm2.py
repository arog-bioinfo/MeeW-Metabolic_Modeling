"""Snakemake local wrapper for CheckM2."""

import os
from pathlib import Path

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
report = Path(snakemake.output.report)
output_directory = report.parent
database_path = Path(
    os.path.expanduser(os.path.expandvars(str(snakemake.params.database_path)))
)
expected_database_file = database_path / "CheckM2_database" / "uniref100.KO.1.dmnd"

protein_extensions = {"faa", "faa.gz", "aa", "aa.gz"}
extension = str(snakemake.params.extension).lstrip(".")
genes_flag = "--genes" if extension in protein_extensions else ""

if database_path.suffix == ".dmnd" or database_path.is_file():
    database_file = database_path
else:
    database_file = expected_database_file

if snakemake.params.download_database and not expected_database_file.exists():
    database_path.mkdir(parents=True, exist_ok=True)
    shell(
        "checkm2 database --download --path {database_path:q} --no_write_json_db {log}"
    )

shell(
    "checkm2 predict "
    " --input {snakemake.input.bins:q} "
    " --output-directory {output_directory:q} "
    " --database_path {database_file:q} "
    " --threads {snakemake.threads} "
    " --extension {extension:q} "
    " {genes_flag} "
    " {snakemake.params.extra} "
    " {log}"
)
