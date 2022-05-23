#########################################
#### creating folder structure #####
#########################################
rule prep:
    input:
        REF=config["reference"],
        BED=config["bed"]

    shell:
        """
        mkdir -p results
        mkdir -p results/PEPPER
        mkdir -p results/bam
        mkdir -p results/report
        mkdir -p results/{SAMPLE}
        mkdir -p results/{SAMPLE}/out
        mkdir -p results/{SAMPLE}/out/QC
        mkdir -p results/{SAMPLE}/out/confident_snv
        mkdir -p results/{SAMPLE}/out/tabular
        mkdir -p results/PEPPER/{SAMPLE}
        mkdir -p results/{SAMPLE}/out/text_file
        mkdir -p results/{output.unique_out_path}
        mkdir -p results/{SAMPLE}/out/snps_to_evaluate/svm_train_data
        mkdir -p results/{SAMPLE}/out/final_vcf_SNIPER
        """

#########################################
#### mapping via minimap2 ###############
#########################################
rule mapping_minimap2:
    input:
        FASTQ=config["data"],
        REF= config["reference"],
        BED= config["bed"]
    output:
        out="results/bam/{SAMPLE}_minimap2.sam",
        bam="results/bam/{SAMPLE}_minimap2.bam",
        bam_sorted="results/bam/{SAMPLE}.fastq_sorted.minimap2.bam",
        filterd_bam="results/bam/{SAMPLE}_filtered.bam"
    threads: workflow.cores * 0.75
    shell:
        """
        minimap2 -B 2 -t {threads} -m 20 -ax map-ont {input.REF} {input.FASTQ}/{SAMPLE}.fastq >  {output.out}
        samtools view -b results/bam/{SAMPLE}_minimap2.sam -o {output.bam}
        samtools sort -o {output.bam_sorted} results/bam/{SAMPLE}_minimap2.bam
        samtools index results/bam/{SAMPLE}.fastq_sorted.minimap2.bam

        #remove unmapped, secondary and other alignment

        samtools view -q 40 -F 2048 -F 256 -b results/bam/{SAMPLE}.fastq_sorted.minimap2.bam > {output.filterd_bam}
        """


