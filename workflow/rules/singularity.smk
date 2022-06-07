########################################################
############ PEPPER called variants#####################
#######################################################

rule PEPPER_calling:
    input:
        bam="results/{SAMPLE}/out/{SAMPLE}_RG_sorted.bam",
        REF=config["reference"],
        
    output:
        OUTPUT_VCF="results/PEPPER/{SAMPLE}/intermediate_files/PEPPER_HP_VARIANT_FULL.vcf.gz"
    #threads: workflow.cores * 0.75
    params:
        flow_cell=config["ont_chemistry"]
    cache: True
    shell:
        """
        singularity exec --bind $PWD envs/pepper_deepvariant_r0.7.sif run_pepper_margin_deepvariant call_variant  --only_pepper -b "{input.bam}" -f "{input.REF}" -o "results/PEPPER/{SAMPLE}" -p "results/PEPPER/{SAMPLE}" -s "{SAMPLE}" -t {threads} --'{params.flow_cell}'
        """

