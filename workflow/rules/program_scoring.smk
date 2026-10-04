rule program_scoring:
    input:
        vst = "results/deseq2/vst_matrix.csv",
        clusters = "results/tables/gene_clusters.csv"
    output:
        scores = "results/tables/program_scores.csv"
    conda:
        "../envs/deseq2.yaml"
    script:
        "../scripts/08_program_scoring.R"