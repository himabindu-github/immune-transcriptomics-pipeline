rule heatmap:
    input:
        vst="results/deseq2/vst_matrix.csv",
        res="results/deseq2/{case}_vs_{control}.csv",
        meta="config/sample_metadata_gsm.tsv"
    output:
        "results/deseq2/{case}_vs_{control}_heatmap.png"
    conda:
        "../envs/deseq2.yaml"
    script:
        "../scripts/05_heatmap.R"
