#!/usr/bin/env Rscript

#require input: svm_training_pepper.txt built database inside the path of args[4]#
set.seed(123)
args = commandArgs(trailingOnly=TRUE)
path=args[1]

#path=args[4]
setwd(path)


library(tidyverse)    # data manipulation and visualization
library(e1071)  
library(tibble)

# specify path to folder containing variant calling files
#INPUT_FOLDER <- "snps_to_evaluate"
OUT_FOLDER <- "pepper_pass"

read_txt_to_df <- function(path){
  read.table(path, sep = '\t', header = TRUE)
}

df <- read_txt_to_df("svm_train_data/svm_training_pepper.txt")

#Preparing for training dataset
df$CHROM = gsub("Pf3D7_","",df$CHROM)
df$CHROM = gsub("_v3","",df$CHROM)
df$REF <- chartr("ACGT", "0123", df$REF)
df$ALT <- chartr("ACGT", "0123", df$ALT)

  
df <-  lapply(df,as.numeric)
df <- as.data.frame(df)
df$Type <- as.factor(df$Type)
str(df)
df<- na.omit(df) 



# tune model to find optimal cost, gamma values
tune.out <- tune(svm, Type~., data = df, kernel = "radial",
                 ranges = list(cost = c(0.1,1,10,100,1000), gamma = c(0.5,1,2,3,4)))

bestmod <- tune.out$best.model
#This is the end of SVM training phase.
#prepare file for SVM prediction, input here is merged file from longshot and pepper unique variants. 
# We will feed them into SVM seperately since each caller has unique error profile.

## returns a df containing the sample IDs and paths of each txt file within the given folder


pred <- function(OUT_FOLDER){
  #setwd(input_path)
  txt_files <- list.files(pattern = "\\.txt$")
  for (file in unique(txt_files)){
    sample_IDs <- gsub("_.*$", "", file)
    file <- read_txt_to_df(file)
    file_pepper <- file[file$ID == 'PEPPER',]
    file_pepper  <- tibble::rownames_to_column(file_pepper, "number")
    file_pepper <- file_pepper[-c(1)]
    file_pepper  <- tibble::rownames_to_column(file_pepper, "index")
    file_pepper_subset <- file_pepper[c("CHROM","POS","REF","ALT","QUAL","AP","GQ","DP","AD","VAF")]
    file_pepper_subset$CHROM <-  gsub("Pf3D7_","",file_pepper_subset$CHROM)
    file_pepper_subset$CHROM <-  gsub("_v3","",file_pepper_subset$CHROM)
    file_pepper_subset$REF <- chartr("ACGT", "0123", file_pepper_subset$REF)
    file_pepper_subset$ALT <- chartr("ACGT", "0123", file_pepper_subset$ALT)
    file_pepper_subset$ALT<- sub('([^,]+\\,).*', '\\1',  file_pepper_subset$ALT)
    file_pepper_subset$ALT<- sub(',', "",  file_pepper_subset$ALT)
    file_pepper_subset$AD<- sub('([^,]+\\,).*', '\\1',  file_pepper_subset$AD)
    file_pepper_subset$AD<- sub(',', "",  file_pepper_subset$AD)
    file_pepper_subset$VAF<- sub('([^,]+\\,).*', '\\1',  file_pepper_subset$VAF)
    file_pepper_subset$VAF<- sub(',', "",  file_pepper_subset$VAF)    
    file_pepper_subset <-  lapply(file_pepper_subset,as.numeric)
    file_pepper_subset <- as.data.frame(file_pepper_subset)
    SVM_filtered <- predict(bestmod, file_pepper_subset)
    SVM_filtered <- as.data.frame(SVM_filtered)
    SVM_filtered$SVM_filtered = gsub("1","PASS",SVM_filtered$SVM_filtered)
    SVM_filtered$SVM_filtered = gsub("2","fail",SVM_filtered$SVM_filtered)
    SVM_filtered <- tibble::rownames_to_column(SVM_filtered, "index")
    out <-  merge(file_pepper,SVM_filtered, by = "index")
    out <- out[out$SVM_filtered == 'PASS',]
    drops <- c("index")
    out <- out[,!(names(out) %in% drops)]
    colnames(out)[which(names(out) == "SVM_filtered")] <- "FILTER"
    colnames(out)[which(names(out) == "CHROM")] <- "#CHROM"
    write.table(out, file = paste0(OUT_FOLDER,"/",sample_IDs,"_pepper_pass.txt"), row.names = FALSE, sep='\t', quote = F)
  }

} 
 
pred(OUT_FOLDER)
