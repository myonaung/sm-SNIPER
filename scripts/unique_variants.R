#!/usr/bin/env Rscript

## Purpose: script to merge tab-delimited txt files
## Description: given a folder of  variant calling txt files (tab delimited),
##              for each sample, we find the set of variants that only appear in
##              one file and not any others ('unique variants'). 
##              Current iteration only handles case where samples have exactly 
##              2 files each.

# specify path to folder containing variant calling files

args = commandArgs(trailingOnly=TRUE)
library(readr)
library(tidyr)
library(tidyverse)    # data manipulation and visualization
library(e1071)  
library(tibble)


INPUT_FOLDER <-  args[1]
OUTPUT_FOLDER <- args[2]

# read tab delimited file into a dataframe
read_txt_to_df <- function(path){
  read.table(path, sep = '\t', header = TRUE)
}

# takes a df containing columns for CHROM and POS, and returns a vector of all
# CHROM POS combinations in the df
get_variants <- function(df){
  paste(df$CHROM, df$POS)
}

# takes in a df containing variant calls (at least columns for GT, ALT and REF)
# and returns a df with the GT column binarised so that any REF alleles are 
# set to 0 and any ALT alleles are set to 1
binarise_genotypes <- function(df){
  # replace any REF alleles in GT with 0
  df$GT <- mapply(function(ref,gt) gsub(ref, "0", gt), df$REF, df$GT)
  # replace any ALT alleles in GT with 1
  df$GT <- mapply(function(alt,gt) gsub(alt, "1", gt), df$ALT, df$GT)
  df
}

# merges file1 and file2, and writes the output to output_path
merge_files <- function(file1_path, file2_path, output_path, sample_id){

  df1 <- read_txt_to_df(file1_path)
  df2 <- read_txt_to_df(file2_path)
  
  # get the variants present in each file
  variants1 <- get_variants(df1)
  variants2 <- get_variants(df2)
  
  # find common variants (which we want to remove)
  variants_to_remove <- intersect(variants1, variants2)
  
  # create new dataframes containing only the unique variants 
  # (can be done in place if required, but left the original dfs untouched in case
  # they need to be available later for some reason)
  df1 <- df1[!(variants1 %in% variants_to_remove),]
  df2 <- df2[!(variants2 %in% variants_to_remove),]
  
  # check if the filtered dfs are empty, since the next section where we add
  # NA columns will break if the dfs are empty
  if (nrow(df1) == 0){
    if (nrow(df2) == 0){
      # if both dfs are empty, nothing to merge so don't write out a file
      # print statement informing which sample wasn't merged
      cat("Merged file for sample", sample_id, 
          "is empty, thus no file will be outputted.", sep=' ')
      return()
    }
    # otherwise if only first df is empty, binarise and write out second df
    df2 <- binarise_genotypes(df2)
    write.table(df2, file = output_path, row.names = FALSE, sep='\t')
    return()
  }
  if (nrow(df2) == 0){
    # if only second df is empty, binarise and write out first df
    df1 <- binarise_genotypes(df1)
    write.table(df1, file = output_path, row.names = FALSE, sep='\t')
    return()
  }
  
  # add columns of NA for any columns that are missing in one df but present in 
  # the other
  df1[setdiff(names(df2), names(df1))] <- NA
  df2[setdiff(names(df1), names(df2))] <- NA
  
  # combine the two dfs into one
  combined_df <- rbind(df1, df2)
  
  # binarise the genotype column
  combined_df <- binarise_genotypes(combined_df)
  
  # write out the combined df to a tab delimited file
  write.table(combined_df, file = output_path, row.names = FALSE, sep='\t', quote = F)
}

# returns a df containing the sample IDs and paths of each txt file within the 
# given folder
get_sample_files <- function(folder_path){
  
  # get a vector of all txt files in the folder
  txt_files <- list.files(folder_path, pattern = "\\.txt$")
  # get a vector of the full paths to the files
  full_paths <- file.path(folder_path, txt_files)
  
  # get a vector of sample ids for each file (replace _ and any text after with
  # empty string)
  sample_IDs <- gsub("_.*$", "", txt_files)
  
  # convert to dataframe
  sample_files <- data.frame(sample_ID = sample_IDs, path = full_paths)
  sample_files
}

# check that the samples have exactly 2 files to merge
# display message to user to inform which files don't meet this criteria, 
# and remove these from the sample_files df. Return a df of only valid samples
validate_sample_files <- function(sample_files){
  
  # find any samples with less or more than 2 associated files
  invalid_samples <- names(which(table(sample_files$sample_ID) != 2))
  
  # if there's at least one invalid file, print message informing user which
  # samples are invalid, and remove these from df
  if(length(invalid_samples)>0){
    cat("The following samples had more or less than 2 associated files and will be ignored:",
        invalid_samples, sep = "\n")
    sample_files <- sample_files[!(sample_files$sample_ID %in% invalid_samples),]
  }
  
  sample_files
}

# get a df of all txt files and corresponding sample names
sample_files <- get_sample_files(INPUT_FOLDER)
# check that the samples are valid and filter out invalid samples
sample_files <- validate_sample_files(sample_files)

# for each unique sample
for (sample in unique(sample_files$sample_ID)){
  # get the files for the sample
  files <- sample_files[(sample_files$sample_ID == sample), 'path']
  merged_file_path <- file.path(OUTPUT_FOLDER, paste0(sample, '_merged.txt'))
  # merge these files. Assumes there is exactly 2 files per sample.
  # This is enforced by the validation step earlier
  merge_files(files[1], files[2], merged_file_path, sample_id = sample)
}
