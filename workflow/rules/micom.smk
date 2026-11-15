rule micom:
    input:
        models=memote_passing_model_files,
        memote_quality=memote_quality_file,
        filtered_models=memote_filtered_models_file,
        abundance=micom_abundance_input,
        medium=lambda wc: selected_medium_output("medium"),
        medium_mapping=lambda wc: selected_medium_output("mapping"),
        medium_status=lambda wc: selected_medium_output("status"),
        scientific_config=simulation_config_inputs,
    output:
        community_growth=f"<results>/micom/{SIMULATION_MEDIUM_TAG}/community_growth.tsv",
        member_growth=f"<results>/micom/{SIMULATION_MEDIUM_TAG}/member_growth.tsv",
        fluxes=f"<results>/micom/{SIMULATION_MEDIUM_TAG}/fluxes.tsv",
        status=f"<results>/micom/{SIMULATION_MEDIUM_TAG}/status.json",
    log:
        f"<logs>/micom/{SIMULATION_MEDIUM_TAG}.log",
    conda:
        "../envs/micom.yaml" if USE_CONDA_ENVS else None
    threads: 1
    params:
        solver=config.get("solver", "cplex"),
        tag=SIMULATION_MEDIUM_TAG,
        mock=config.get("execution", {}).get("mock", False),
    script:
        "../scripts/micom_script.py"
