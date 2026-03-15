rule fastqc:
    input:
        "data/raw/{sample}.fastq"
    output:
        "results/qc/{sample}_fastqc.html"
    conda:
        "../envs/fastqc.yaml"
    shell:
        """
        fastqc {input} -o results/qc
        """
