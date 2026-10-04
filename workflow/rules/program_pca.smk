rule program_pca:
    input:
        scores = "results/tables/program_scores.csv",
        meta   = "config/sample_metadata_gsm.tsv"
    output:
        ppca = "results/figures/gsva/program_pca.png"
    conda:
        "../envs/deseq2.yaml"
    script:
        "../scripts/09_program_pca.R"