
rule deseq2_global:
    input:
        counts="data/geo_counts/cd8_counts.tsv",
        meta="config/sample_metadata_gsm.tsv"
    output:
        vst="results/deseq2/vst_matrix.csv",
        pca="results/deseq2/global_PCA_plot.pdf",
        dds="results/deseq2/dds.rds"
    conda:
        "../envs/deseq2.yaml"
    script:
        "../scripts/02-global-deseq2.R"
