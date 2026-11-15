rule discover_individual_medium:
    input:
        model="<results>/gems/{mag}.xml",
        filtered_models=memote_filtered_models_file,
        scientific_config=simulation_config_inputs,
    output:
        medium="<results>/solo/{mag}/minimal_medium.csv",
        summary="<results>/solo/{mag}/summary.json",
    log:
        "<logs>/solo/{mag}.log",
    conda:
        "../envs/micom.yaml" if USE_CONDA_ENVS else None
    threads: 1
    params:
        solver=config.get("solver", "cplex"),
        min_growth=simulation_medium_config["min_growth"],
        open_exchanges=simulation_medium_config["open_exchanges"],
        exchange_bound=simulation_medium_config["exchange_bound"],
        mock=config.get("execution", {}).get("mock", False),
    script:
        "../scripts/discover_medium.py"


rule prepare_simulation_medium:
    input:
        models=memote_passing_model_files,
        filtered_models=memote_filtered_models_file,
        recipes=selected_recipe_inputs,
        predefined=selected_predefined_medium_input,
        scientific_config=simulation_config_inputs,
    output:
        medium=f"<results>/media/{SIMULATION_MEDIUM_TAG}.csv",
        mapping=f"<results>/media/{SIMULATION_MEDIUM_TAG}.mapping.tsv",
        status=f"<results>/media/{SIMULATION_MEDIUM_TAG}.status.json",
    log:
        f"<logs>/media/{SIMULATION_MEDIUM_TAG}.log",
    conda:
        "../envs/micom.yaml" if USE_CONDA_ENVS else None
    threads: 1
    params:
        source=simulation_medium_config["source"],
        aggregation=simulation_medium_config["aggregation"],
        predefined_label=simulation_predefined_label,
        tag=SIMULATION_MEDIUM_TAG,
        min_growth=simulation_medium_config["min_growth"],
        open_exchanges=simulation_medium_config["open_exchanges"],
        exchange_bound=simulation_medium_config["exchange_bound"],
        mock=config.get("execution", {}).get("mock", False),
    script:
        "../scripts/prepare_simulation_medium.py"
