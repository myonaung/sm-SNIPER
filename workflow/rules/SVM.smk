#################################
#### Support Vector Machine #####
#################################
rule SVM:
    input:
        unique_out="results/{SAMPLE}/out/snps_to_evaluate/{SAMPLE}_merged.txt",
        path_SVM="results/{SAMPLE}/out/snps_to_evaluate"
    output:
        longshot_pass="results/{SAMPLE}/out/snps_to_evaluate/longshot_pass/{SAMPLE}_longshot_pass.txt",
        pepper_pass="results/{SAMPLE}/out/snps_to_evaluate/pepper_pass/{SAMPLE}_pepper_pass.txt"
    shell:
        """
        scp resources/svm_training_longshot.txt results/{SAMPLE}/out/snps_to_evaluate/svm_train_data
        scp resources/svm_training_pepper.txt results/{SAMPLE}/out/snps_to_evaluate/svm_train_data
        scripts/svm_longshot.R results/{SAMPLE}/out/snps_to_evaluate/
        scripts/svm_pepper.R results/{SAMPLE}/out/snps_to_evaluate/
        """