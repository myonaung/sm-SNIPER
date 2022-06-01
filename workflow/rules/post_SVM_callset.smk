#################################
#### Post SVM Processing ########
#################################
rule post_SVM:
    input:
        longshot_pass="results/{SAMPLE}/out/snps_to_evaluate/longshot_pass/{SAMPLE}_longshot_pass.txt",
        pepper_pass="results/{SAMPLE}/out/snps_to_evaluate/pepper_pass/{SAMPLE}_pepper_pass.txt"
    output:
        longshot_pass_sites="results/{SAMPLE}/out/snps_to_evaluate/longshot_pass/{SAMPLE}_longshot_svm_pass_sites",
        pepper_pass_sites="results/{SAMPLE}/out/snps_to_evaluate/pepper_pass/{SAMPLE}_pepper_svm_pass_sites",
    shell:
        """
        scripts/post_SVM.sh {input.longshot_pass} > {output.longshot_pass_sites}
        scripts/post_SVM.sh {input.pepper_pass} > {output.pepper_pass_sites}
        """

rule svm_longshot_pass:
    input:
        longshot_amplicon= "results/{SAMPLE}/out/{SAMPLE}_amplicon.vcf.gz",
        longshot_pass_sites="results/{SAMPLE}/out/snps_to_evaluate/longshot_pass/{SAMPLE}_longshot_svm_pass_sites"
    output:
        longshot_SVM_pass="results/{SAMPLE}/out/{SAMPLE}_longshot_SVM_pass.vcf"
    shell:
        """
        vcftools --gzvcf {input.longshot_amplicon} --positions {input.longshot_pass_sites} --recode --out results/{SAMPLE}/out/longshot_final
        mv results/{SAMPLE}/out/longshot_final.recode.vcf {output.longshot_SVM_pass}
        """

rule final_longshot:
    input:
        longshot_SVM_pass="results/{SAMPLE}/out/{SAMPLE}_longshot_SVM_pass.vcf",
        confident_snv_index="results/{SAMPLE}/out/confident_snv/{SAMPLE}_confident_snps.vcf.gz"

    output:
        index_longshot_SVM_pass="results/{SAMPLE}/out/{SAMPLE}_longshot_SVM_pass.vcf.gz",
        final_longshot="results/{SAMPLE}/out/{SAMPLE}_final_longshot.vcf",
        final_longshot_index="results/{SAMPLE}/out/{SAMPLE}_final_longshot.vcf.gz"
    shell:
        """
        bgzip -c {input.longshot_SVM_pass} > {output.index_longshot_SVM_pass}
        tabix {output.index_longshot_SVM_pass}

        bcftools concat -a results/{SAMPLE}/out/{SAMPLE}_longshot_SVM_pass.vcf.gz {input.confident_snv_index} -o {output.final_longshot}

        bgzip -c results/{SAMPLE}/out/{SAMPLE}_final_longshot.vcf > {output.final_longshot_index}
        tabix {output.final_longshot_index}
        """

rule svm_PEPPER:
    input:
        index_pepper_amplicon= "results/{SAMPLE}/out/{SAMPLE}_PEPPER_VARIANT_FULL.vcf.gz",
        pepper_pass_sites="results/{SAMPLE}/out/snps_to_evaluate/pepper_pass/{SAMPLE}_pepper_svm_pass_sites"
    output:
        pepper_SVM_pass="results/{SAMPLE}/out/{SAMPLE}_PEPPER_SVM_pass.vcf",
        pepper_SVM_pass_index="results/{SAMPLE}/out/{SAMPLE}_PEPPER_SVM_pass.vcf.gz"
    shell:
        """
        vcftools --gzvcf {input.index_pepper_amplicon} --positions {input.pepper_pass_sites} --recode --out results/{SAMPLE}/out/pepper_final
        mv results/{SAMPLE}/out/pepper_final.recode.vcf {output.pepper_SVM_pass}

        bgzip -c results/{SAMPLE}/out/{SAMPLE}_PEPPER_SVM_pass.vcf > {output.pepper_SVM_pass_index}
        tabix {output.pepper_SVM_pass_index}
        """

rule format_conversion:
    input:
        bam="results/bam/{SAMPLE}.fastq_sorted.minimap2.bam",
        REF= config["reference"],
        pepper_SVM_pass_index="results/{SAMPLE}/out/{SAMPLE}_PEPPER_SVM_pass.vcf.gz"
    output:
        pepper_to_longshot="results/{SAMPLE}/out/{SAMPLE}_pepper_to_longshot.vcf",
        pepper_to_longshot_index="results/{SAMPLE}/out/{SAMPLE}_pepper_to_longshot.vcf.gz"
    shell:
        """
        longshot -b {input.bam} -f {input.REF} -v {input.pepper_SVM_pass_index} -o {output.pepper_to_longshot} --sample_id {SAMPLE} -E 0.01 -e 2 -C 5000 -c 1 -F -P 0.00001 -n --output-ref  -q 5 -B 50 -l 20 -y 10 -a 10
        
        bgzip -c results/{SAMPLE}/out/{SAMPLE}_pepper_to_longshot.vcf > {output.pepper_to_longshot_index}
        tabix {output.pepper_to_longshot_index}
        """
##################################
#### creating final variants #####
##################################

rule final_variants:
    input:
        final_longshot_index="results/{SAMPLE}/out/{SAMPLE}_final_longshot.vcf.gz",
        pepper_to_longshot_index="results/{SAMPLE}/out/{SAMPLE}_pepper_to_longshot.vcf.gz"
    output:
        final_variants="results/{SAMPLE}/out/final_vcf_SNIPER/{SAMPLE}_final_SNIPER.vcf"
    shell:
        """
        bcftools concat -a {input.final_longshot_index} {input.pepper_to_longshot_index} -o {output.final_variants}
        """

