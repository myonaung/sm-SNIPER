#!/usr/bin/env Rscript

#require input: svm_training_pepper.txt built database inside the path of args[4]#

set.seed(123)
args = commandArgs(trailingOnly=TRUE)
path=args[1]


setwd(path)

library(tidyverse)    # data manipulation and visualization
library(e1071)  
library(tibble)

# specify path to folder containing variant calling files
#INPUT_FOLDER <- "out/snps_to_evaluate"
OUT_FOLDER <- "longshot_pass"

read_txt_to_df <- function(path){
  read.table(path, sep = '\t', header = TRUE)
}

df <- read_txt_to_df("svm_train_data/svm_training_longshot.txt")

#Preparing for training dataset
df$CHROM = gsub("Pf3D7_","",df$CHROM)
df$CHROM = gsub("_v3","",df$CHROM)
df$REF <- chartr("ACGT", "0123", df$REF)
df$ALT <- chartr("ACGT", "0123", df$ALT)
df$SC <- chartr("ACGT", "0123", df$SC)


df <-  lapply(df,as.numeric)
df <- as.data.frame(df)
df$Type <- as.factor(df$Type)
str(df)
df<- na.omit(df) 

# tune model to find optimal cost, gamma values
tune.out <- tune(svm, Type~., data = df, kernel = "radial",
                 ranges = list(cost = c(0.1,1,10,100,1000), gamma = c(0.5,1,2,3,4)))

#This is the end of SVM training phase.
#prepare file for SVM prediction, input here is merged file from longshot and pepper unique variants. 
# We will feed them into SVM seperately since each caller has unique error profile.

## returns a df containing the sample IDs and paths of each txt file within the given folder
bestmod <- tune.out$best.model

pred <- function(OUT_FOLDER){
  #setwd(input_path)
  txt_files <- list.files(pattern = "\\.txt$")
  for (file in unique(txt_files)){
    sample_IDs <- gsub("_.*$", "", file)
    file <- read_txt_to_df(file)
    file_longshot <- file[file$ID == 'longshot',]
    if(nrow(file_longshot) == 0) return(NULL)
    file_longshot  <- tibble::rownames_to_column(file_longshot, "number")
    file_longshot <- file_longshot[-c(1)]
    file_longshot  <- tibble::rownames_to_column(file_longshot, "index")
    file_longshot_subset <- file_longshot[c("CHROM","POS","REF","ALT","QUAL","SC","AQ","DP","GQ","UQ")]
    file_longshot_subset$CHROM <-  gsub("Pf3D7_","",file_longshot_subset$CHROM)
    file_longshot_subset$CHROM <-  gsub("_v3","",file_longshot_subset$CHROM)
    file_longshot_subset$REF <- chartr("ACGT", "0123", file_longshot_subset$REF)
    file_longshot_subset$ALT <- chartr("ACGT", "0123", file_longshot_subset$ALT)
   file_longshot_subset$SC <- chartr("ACGT", "0123", file_longshot_subset$SC) 
   file_longshot_subset <-  lapply(file_longshot_subset,as.numeric)
    file_longshot_subset <- as.data.frame(file_longshot_subset)
    SVM_filtered <- predict(bestmod, file_longshot_subset)
    SVM_filtered <- as.data.frame(SVM_filtered)
    SVM_filtered$SVM_filtered = gsub("1","PASS",SVM_filtered$SVM_filtered)
    SVM_filtered$SVM_filtered = gsub("2","fail",SVM_filtered$SVM_filtered)
    SVM_filtered <- tibble::rownames_to_column(SVM_filtered, "index")
    out <-  merge(file_longshot,SVM_filtered, by = "index")
    out <- out[out$SVM_filtered == 'PASS',]
    drops <- c("index")
    out <- out[,!(names(out) %in% drops)]
    colnames(out)[which(names(out) == "SVM_filtered")] <- "FILTER"
    colnames(out)[which(names(out) == "CHROM")] <- "#CHROM"
    write.table(out, file = paste0(OUT_FOLDER,"/",sample_IDs,"_longshot_pass.txt"), row.names = FALSE, sep='\t', quote = F)
  }
}


pred(OUT_FOLDER)
