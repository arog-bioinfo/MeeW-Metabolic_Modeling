rule stage_checkm2_faa:
    input:
        faas=checkm2_faa_inputs,
        direct_inputs_helper=direct_inputs_helper_input,
    output:
        bins=directory("<results>/qa/checkm2/input"),
    log:
        "<logs>/qa/checkm2/stage_faa.log",
    conda:
        "../envs/python.yaml"
    params:
        faa_by_mag=checkm2_faa_by_mag_param,
        extension=checkm2_config.get("extension", "faa"),
        bakta_dir=standalone_bakta_dir_param,
        direct_annotation=direct_annotation_param,
    script:
        "../scripts/stage_checkm2_faa.py"


rule checkm2:
    input:
        bins="<results>/qa/checkm2/input",
    output:
        report=checkm2_report,
    log:
        "<logs>/qa/checkm2/checkm2.log",
    conda:
        "../envs/checkm2.yaml" if USE_CONDA_ENVS else None
    threads: checkm2_config.get("threads", 8)
    params:
        database_path=checkm2_config.get("database_path", ""),
        download_database=checkm2_config.get("download_database", False),
        extension=checkm2_config.get("extension", "faa"),
        extra=checkm2_config.get("extra", ""),
    script:
        "../scripts/checkm2.py"


checkpoint external_qa_filter:
    input:
        faa_samples=faa_sample_sheet_input,
        bakta_done=bakta_done_input,
        annotation_dirs=protein_source_directories_input,
        qa_tsv=(
            checkm2_report
            if qa_uses_checkm2
            else qa_tsv if qa_enabled and qa_tsv else []
        ),
    output:
        summary="<results>/qa/external_qa_summary.tsv",
        filtered="<results>/qa/filtered_mags.txt",
    log:
        "<logs>/qa/external_qa_filter.log",
    conda:
        "../envs/python.yaml"
    params:
        faa_by_mag=external_qa_faa_by_mag_param,
        mags=external_qa_mags_param,
        bakta_dir=standalone_bakta_dir_param,
        direct_annotation=direct_annotation_param,
        enabled=qa_enabled,
        id_col=qa_id_col,
        completeness_col=qa_completeness_col,
        contamination_col=qa_contamination_col,
        min_completeness=qa_min_completeness,
        max_contamination=qa_max_contamination,
    script:
        "../scripts/filter_external_qa.py"
