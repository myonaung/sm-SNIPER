#########################################
#### haplotype calling  ################
#########################################
rule phase:
    input:
        final_variants="results/{SAMPLE}/out/final_vcf_SNIPER/{SAMPLE}_final_SNIPER.vcf",
        REF= config["reference"],
        filterd_bam="results/bam/{SAMPLE}_filtered.bam"
    output:
        out_phase="results/{SAMPLE}/out/final_vcf_SNIPER/{SAMPLE}_final_SNIPER_phased.vcf",
    shell:
        """
        whatshap phase -o {output.out_phase} --reference {input.REF} {input.final_variants} {input.filterd_bam} --ignore-read-groups --mapq 40 
        """


