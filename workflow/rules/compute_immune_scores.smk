rule compute_immune_scores:
    input:
        vst="results/deseq2/vst_matrix.csv",
        signatures="resources/gene_signatures.csv"
    output:
        scores="results/tables/immune_scores.csv"
    script:
        "../scripts/13b_immune_scores.R"