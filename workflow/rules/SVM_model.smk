#########################################
#### SVM model generation #####
#########################################
rule SVM_model:
    input:
        longshot_truth="resources/svm_training_longshot.txt",
        pepper_truth="resources/svm_training_pepper.txt",

    output:
        longshot_model="resources/longshotmodel.rds",
        pepper_model="resources/peppermodel.rds"

    shell:
        """
        scripts/svm_longshot_model.R resources/svm_training_longshot.txt
        scripts/svm_pepper_model.R resources/svm_training_pepper.txt
        """
