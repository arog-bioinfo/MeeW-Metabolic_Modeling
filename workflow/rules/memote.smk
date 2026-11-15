rule memote:
    input:
        model="<results>/gems/{mag}.xml",
    output:
        report="<results>/memote/{mag}/report.html",
        score="<results>/memote/{mag}/score.json",
    log:
        "<logs>/memote/{mag}.log",
    conda:
        "../envs/memote.yaml" if USE_CONDA_ENVS else None
    params:
        mock=config.get("execution", {}).get("mock", False),
    script:
        "../scripts/memote.py"


rule memote_summary:
    input:
        scores=first_qa_passing_memote_scores,
    output:
        summary="<results>/memote/memote_summary.tsv",
    log:
        "<logs>/memote/memote_summary.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/memote_summary.py"


checkpoint filter_memote:
    input:
        scores=first_qa_passing_memote_scores,
    output:
        quality="<results>/memote/memote_quality.tsv",
        filtered="<results>/memote/filtered_models.txt",
    log:
        "<logs>/memote/filter_memote.log",
    conda:
        "../envs/python.yaml"
    params:
        min_total_score=config.get("memote", {}).get("min_total_score", 0.4),
    script:
        "../scripts/filter_memote.py"
