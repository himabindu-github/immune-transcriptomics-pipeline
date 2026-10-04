rule qc_contamination_check:
    input:
        vst = "results/deseq2/vst_matrix.csv",
        meta = "config/sample_metadata_gsm.tsv"
    output:
        table = "results/qc/contamination_scores.csv",
        plot = "results/qc/neutrophil_contamination_check.png"
    conda:
        "../envs/deseq2.yaml"
    script:
        "../scripts/qc_contamination_check.R"
