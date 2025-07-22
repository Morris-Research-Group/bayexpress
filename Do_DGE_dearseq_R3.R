library(tidyverse)
library(dearseq)
library(edgeR)
library(data.table)
library(matrixStats)

args = commandArgs(trailingOnly=TRUE)

do_DGE <- function(input_file, output_file_name){
  
  input_data <- read.csv(input_file, header = TRUE, sep = ",")
  count_data <- input_data %>% dplyr::select(-locus_name)

  group <- data.matrix(c(1,1,1,2,2,2))
  count_data <- as.matrix(count_data)

  res_genes <- dearseq::dear_seq(exprmat = count_data,
                                 variables2test=group,
                                 which_test='permutation',
                                 preprocessed = FALSE, n_perm=700)
  
  write.csv(as.data.frame(res_genes$pvals), 
            file=output_file_name)
}

do_DGE(args[1], args[2])