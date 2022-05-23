#########################################
#### quality check to alignment files #####
#########################################
rule QC:
    input:
        FASTQ=config["data"],
        bam="results/bam/{SAMPLE}.fastq_sorted.minimap2.bam",
        REF=config["reference"],
        BED=config["bed"]
    output:
        total_depth="results/{SAMPLE}/out/QC/{SAMPLE}_amplicon_total_depth.txt",
        amplicon_depth="results/{SAMPLE}/out/QC/{SAMPLE}_amplicon_depth.txt",
        mapping_summary="results/{SAMPLE}/out/QC/{SAMPLE}_mapping_summary.txt",
        sample_mean_depth="results/{SAMPLE}/out/QC/{SAMPLE}_mean_depth.txt",

    shell:
        """
        samtools bedcov {input.BED} {input.bam} > {output.total_depth}
        samtools depth -b {input.BED} {input.bam} > {output.amplicon_depth}
        samtools flagstats {input.bam} > {output.mapping_summary}

        awk -f scripts/average_depth.awk results/{SAMPLE}/out/QC/{SAMPLE}_amplicon_depth.txt | cut -d" " -f3> {output.sample_mean_depth}
        """
