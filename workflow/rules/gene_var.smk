rule gene_variability:
    input:
        vst = "results/deseq2/vst_matrix.csv",   
        meta = "config/sample_metadata_gsm.tsv"  
    output:
        var_genes = "results/tables/genes_var.txt"

    conda:
        "../envs/deseq2.yaml"

    script:
        "../scripts/gene_variability.R"
