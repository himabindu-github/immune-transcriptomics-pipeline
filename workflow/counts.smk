rule download_gtf:
    output:
        "reference/annotation.gtf.gz"
    shell:
        """
        wget -O {output} \
        https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz
        """


rule featurecounts:
    input:
        bam="results/alignment/{sample}.bam",
        gtf="reference/annotation.gtf.gz"
    output:
        "results/counts/{sample}_counts.txt"
    threads: 2
    conda:
        "../envs/featurecounts.yaml"
    shell:
        """
        mkdir -p results/counts

        featureCounts \
            -T {threads} \
            -a {input.gtf} \
            -o {output} \
            {input.bam}
        """
