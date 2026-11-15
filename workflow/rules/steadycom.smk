rule steadycom:
    input:
        models=memote_passing_model_files,
        filtered_models=memote_filtered_models_file,
        medium=lambda wc: selected_medium_output("medium"),
        medium_mapping=lambda wc: selected_medium_output("mapping"),
        medium_status=lambda wc: selected_medium_output("status"),
        scientific_config=simulation_config_inputs,
    output:
        community_growth=f"<results>/steadycom/{SIMULATION_MEDIUM_TAG}/community_growth.tsv",
        abundances=f"<results>/steadycom/{SIMULATION_MEDIUM_TAG}/abundances.tsv",
        exchanges=f"<results>/steadycom/{SIMULATION_MEDIUM_TAG}/exchanges.tsv",
        variability=f"<results>/steadycom/{SIMULATION_MEDIUM_TAG}/variability.tsv",
        status=f"<results>/steadycom/{SIMULATION_MEDIUM_TAG}/status.json",
    log:
        f"<logs>/steadycom/{SIMULATION_MEDIUM_TAG}.log",
    conda:
        "../envs/steadycom.yaml" if USE_CONDA_ENVS else None
    threads: 1
    params:
        solver=config.get("solver", "cplex"),
        objective=config.get("steadycom", {}).get("objective", "community_growth"),
        tag=SIMULATION_MEDIUM_TAG,
        variability_fraction=0.9,
        mock=config.get("execution", {}).get("mock", False),
    script:
        "../scripts/steadycom.py"
