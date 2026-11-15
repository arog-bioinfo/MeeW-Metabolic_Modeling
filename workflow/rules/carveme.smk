rule initialize_carveme_index:
    input:
        helper=_CARVEME_INDEX_HELPER_SOURCE,
    output:
        ready=CARVEME_INDEX_READINESS,
    log:
        CARVEME_INDEX_LOG,
    conda:
        "../envs/carveme.yaml"
    script:
        "../scripts/initialize_carveme_index.py"


rule carveme:
    input:
        genome=faa_for_mag,
        qa="<results>/qa/filtered_mags.txt",
        faa_samples=faa_sample_sheet_input,
        gtdb_classification=gtdb_classification_input,
        index_ready=carveme_index_readiness_input,
    output:
        model="<results>/gems/{mag}.xml",
    log:
        "<logs>/carveme/{mag}.log",
    conda:
        "../envs/carveme.yaml" if USE_CONDA_ENVS else None
    threads: 4
    params:
        solver=config.get("solver", "scip"),
        media=config.get("media", "M9"),
        mock=config.get("execution", {}).get("mock", False),
        direct_annotation=direct_annotation_param,
    script:
        "../scripts/carveme.py"
