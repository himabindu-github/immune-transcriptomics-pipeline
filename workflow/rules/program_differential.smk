rule program_differential:
    input:
        "results/tables/program_scores.csv",
        "config/sample_metadata_gsm.tsv"
    output:
        "results/tables/program_differential.csv"
    conda:
        "../envs/deseq2.yaml"
    script:
        "../scripts/10_program_differential.R"