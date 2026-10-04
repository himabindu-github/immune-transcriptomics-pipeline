rule disease_signatures:
    input:
        summary="results/tables/program_summary.csv",
        annotations="results/tables/program_annotations.csv"
    output:
        signatures="results/tables/disease_signatures.csv"
    script:
        "../scripts/12_disease_signatures.R"