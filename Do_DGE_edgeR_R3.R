#!/usr/bin/env Rscript
library(tidyverse)
library(edgeR)
library(limma)
library(data.table)


args = commandArgs(trailingOnly=TRUE)


do_DGE <- function(input_file, output_file_name, min_cpm=0.0){

  input_data <- read.csv(input_file, header = TRUE, sep = ",")
  count_data <- input_data %>% dplyr::select(-locus_name)
  
  # creating group array? factor? thingy
  group <- factor(c(1,1,1,2,2,2))
  
  # creating edgeR DGEList object
  dgel <- edgeR::DGEList(
    counts = count_data,
    group = group,
    genes = input_data$locus_name
  )
  
  # filtering for min_cpm
  print(paste("Row count before CPM filtering =", nrow(dgel)))
  
  keep <- rowSums(edgeR::cpm(dgel) >= min_cpm) == length(group)
  dgel <- dgel[keep, , keep.lib.sizes = FALSE]
  
  print(paste("Row count before CPM filtering =", nrow(dgel)))
  
  # Also filtering using edgeR's filtering function
  keep <- filterByExpr(dgel)
  dgel <- dgel[keep,,keep.lib.sizes=FALSE]
  
  # calculating library sizes
  dgel <- calcNormFactors(dgel)
  
  # idk why this is necessary but I am creating design matrix anyway
  design <- model.matrix(~group)
  
  # estimating dispersion is required prior to glm fit
  dgel <- estimateDisp(dgel,design, robust = TRUE)
  
  # finally doing the DGE 
  fit <- glmQLFit(dgel,design)
  qlf <- glmQLFTest(fit,coef=ncol(fit$design))
  
  # writing the final output to csv file in results
  top <- edgeR::topTags(qlf, n = Inf)
  
  write.csv(as.data.frame(top), file=output_file_name, quote = FALSE, row.names = FALSE)
}

do_DGE(args[1], args[2], args[3])





