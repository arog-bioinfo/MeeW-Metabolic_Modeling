rule smetana_global:
    input:
        models=memote_passing_model_files,
        filtered_models=memote_filtered_models_file,
        media_db=smetana_media_db_input,
    output:
        global_scores="<results>/smetana/global.tsv",
    log:
        "<logs>/smetana/global.log",
    benchmark:
        "<logs>/smetana/global.benchmark.tsv"
    conda:
        "../envs/smetana.yaml" if USE_CONDA_ENVS else None
    params:
        mode="global",
        media=smetana_config.get("media", []),
        solver=config.get("solver", "scip"),
        mock=config.get("execution", {}).get("mock", False),
    script:
        "../scripts/smetana.py"


rule smetana_detailed:
    input:
        models=memote_passing_model_files,
        filtered_models=memote_filtered_models_file,
        media_db=smetana_media_db_input,
    output:
        detailed="<results>/smetana/detailed.tsv",
    log:
        "<logs>/smetana/detailed.log",
    benchmark:
        "<logs>/smetana/detailed.benchmark.tsv"
    conda:
        "../envs/smetana.yaml" if USE_CONDA_ENVS else None
    params:
        mode="detailed",
        media=smetana_config.get("media", []),
        solver=config.get("solver", "scip"),
        mock=config.get("execution", {}).get("mock", False),
    script:
        "../scripts/smetana.py"
