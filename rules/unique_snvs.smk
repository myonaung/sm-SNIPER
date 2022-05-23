###################################################
#### evaluating single caller supported snvs #####
###################################################
rule feature_extraction:
    input:
        longshot_table_f="results/{SAMPLE}/out/text_file/{SAMPLE}_longshot.txt",
        PEPPER_table_f="results/{SAMPLE}/out/text_file/{SAMPLE}_PEPPER.txt"
    output:
        unique_out="results/{SAMPLE}/out/snps_to_evaluate/{SAMPLE}_merged.txt"
    shell:
        """
        scripts/unique_variants.R results/{SAMPLE}/out/text_file results/{SAMPLE}/out/snps_to_evaluate
        """