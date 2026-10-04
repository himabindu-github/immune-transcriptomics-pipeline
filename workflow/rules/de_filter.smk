rule de_filtering:
    input:
        var_genes = "results/tables/genes_var.txt",
        res = "results/deseq2/{case}_vs_{control}.csv"
    output:
        var_de_genes = "results/tables/{case}_vs_{control}_genes_var_de.txt",
        de_genes = "results/tables/{case}_vs_{control}_de_genes.txt"
    conda:
        "../envs/deseq2.yaml"
    script:
        "../scripts/de_filtering.R"
