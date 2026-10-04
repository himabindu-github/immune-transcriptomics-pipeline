rule deseq2_compare:
    input:
        dds="results/deseq2/dds.rds"

    output:
        "results/deseq2/{case}_vs_{control}.csv"

    params:
        case=lambda wc: wc.case,
        control=lambda wc: wc.control

    conda:
        "../envs/deseq2.yaml"

    script:
        "../scripts/02-deseq2-contrast.R"
