rule run_gsva:
    input:
        vst="results/deseq2/vst_matrix.csv",
        signatures="resources/gene_signatures.csv"
    output:
        gsva_scores="results/tables/gsva_scores.csv"
    
    script:
        "../scripts/compute_gsva_scores.R"