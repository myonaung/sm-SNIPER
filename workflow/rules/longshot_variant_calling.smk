#########################################
#### Longshot variant calling step #####
#########################################
rule longshot_calling:
    input:
        bam="results/bam/{SAMPLE}.fastq_sorted.minimap2.bam",
        REF= config["reference"],
        BED= config["bed"]
    output:
        longshot="results/{SAMPLE}/out/{SAMPLE}.vcf",
        longshot_amplicon="results/{SAMPLE}/out/{SAMPLE}_amplicon.vcf"
    params:
        min_c=config["min_coverage"],
        min_alt_frac=config["min_alt_frac"],
        max_c=config["max_coverage"]
    threads: workflow.cores * 0.75
    shell:
        """
        longshot -S -b {input.bam} -f {input.REF} -o {output.longshot} -c {params.min_c} -F -E {params.min_alt_frac} -C {params.max_c} -I 3 --sample_id {SAMPLE}
        intersectBed -a {output.longshot} -b {input.BED} -header > {output.longshot_amplicon}
        """

