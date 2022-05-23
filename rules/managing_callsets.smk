rule post_PEPPER:
    input:
        BED=config["bed"]
    output:
        PEPPER_amplicon="results/{SAMPLE}/out/{SAMPLE}_PEPPER_VARIANT_FULL.vcf"
    shell:
        """
        intersectBed -a results/PEPPER/{SAMPLE}/PEPPER_VARIANT_FULL.vcf.gz -b {input.BED} -header > {output.PEPPER_amplicon}
        """

#########################################
#### Building confident callset ########
#########################################

rule confident_snvs:
    input:
        longshot_amplicon="results/{SAMPLE}/out/{SAMPLE}_amplicon.vcf",
        PEPPER_amplicon="results/{SAMPLE}/out/{SAMPLE}_PEPPER_VARIANT_FULL.vcf"
    output:
        index_longshot_amplicon= "results/{SAMPLE}/out/{SAMPLE}_amplicon.vcf.gz",
        index_pepper_amplicon="results/{SAMPLE}/out/{SAMPLE}_PEPPER_VARIANT_FULL.vcf.gz",
        confident_snv_index="results/{SAMPLE}/out/confident_snv/{SAMPLE}_confident_snps.vcf.gz"
    shell:
        """
        bgzip -c {input.longshot_amplicon} > {output.index_longshot_amplicon}
        tabix  {output.index_longshot_amplicon}
        bgzip -c {input.PEPPER_amplicon} > {output.index_pepper_amplicon}
        tabix   {output.index_pepper_amplicon}
        
        bcftools isec -c all results/{SAMPLE}/out/{SAMPLE}_amplicon.vcf.gz results/{SAMPLE}/out/{SAMPLE}_PEPPER_VARIANT_FULL.vcf.gz -p results/{SAMPLE}/out/confident_snv/{SAMPLE} -n=2 -w 1

        bgzip -c results/{SAMPLE}/out/confident_snv/{SAMPLE}/0000.vcf > {output.confident_snv_index}
        tabix {output.confident_snv_index}
        """

###############################################################
#### converting vcf to table format for downstream SVM ########
##############################################################

rule tabular:
    input:
        longshot_amplicon="results/{SAMPLE}/out/{SAMPLE}_amplicon.vcf",
        PEPPER_amplicon="results/{SAMPLE}/out/{SAMPLE}_PEPPER_VARIANT_FULL.vcf"
    output:
        longshot_table="results/{SAMPLE}/out/tabular/{SAMPLE}_longshot_amplicon.txt",
        PEPPER_table="results/{SAMPLE}/out/tabular/{SAMPLE}_PEPPER_amplicon.txt"
    shell:
        """
        gatk VariantsToTable -V {input.longshot_amplicon} -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F AC -F AQ -F SC -F PH -GF GT -F DP -GF GQ -GF UQ -O {output.longshot_table}  --show-filtered
        gatk VariantsToTable -V {input.PEPPER_amplicon} -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -GF GT -GF AP -GF GQ -GF DP -GF AD -GF VAF -O {output.PEPPER_table}  --show-filtered
        """

rule file_modify:
    input:
        input_1="results/{SAMPLE}/out/tabular/{SAMPLE}_longshot_amplicon.txt",
        input_2="results/{SAMPLE}/out/tabular/{SAMPLE}_PEPPER_amplicon.txt"
    output: 
        longshot_table_f="results/{SAMPLE}/out/text_file/{SAMPLE}_longshot.txt",
        PEPPER_table_f="results/{SAMPLE}/out/text_file/{SAMPLE}_PEPPER.txt"
    shell:
        """
        scripts/longshot_tabular.sh {input.input_1} {SAMPLE} > {output.longshot_table_f}
        scripts/pepper_tabular.sh {input.input_2} {SAMPLE} > {output.PEPPER_table_f}
        """        


