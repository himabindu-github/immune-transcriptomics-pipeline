rule align_reads:
    input:
        r1="data/raw/{sample}_1.fastq.gz",
        r2="data/raw/{sample}_2.fastq.gz",
    output:
        bam="results/alignment/{sample}.bam"
    threads: 2
    params:
        index="reference/genome_index/genome"
    conda:
        "../envs/hisat2.yaml"
    shell:
        """
        hisat2 -x {params.index} \
               -1 {input.r1} \
               -2 {input.r2} \
               -p {threads} \
        | samtools view -bS - \
        | samtools sort -@ {threads} -o {output.bam}
        
        """


rule index_bam:
    input:
        "results/alignment/{sample}.bam"
    output:
        "results/alignment/{sample}.bam.bai"
    conda:
        "../envs/hisat2.yaml"
    shell:
        """
        samtools index {input}
        """


rule flagstat:
    input:
        "results/alignment/{sample}.bam"
    output:
        "results/qc/{sample}_flagstat.txt"
    conda:
        "../envs/hisat2.yaml"
    shell:
        """
        samtools flagstat {input} > {output}
        """
