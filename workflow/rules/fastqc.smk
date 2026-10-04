rule fastqc:
    input:
        r1="data/raw/{sample}_1.fastq.gz",
        r2="data/raw/{sample}_2.fastq.gz"
    output:
        "results/qc/{sample}_1_fastqc.html",
        "results/qc/{sample}_2_fastqc.html"
    conda:
        "../envs/fastqc.yaml"
    shell:
        """
        fastqc {input.r1} {input.r2} -o results/qc
        """
