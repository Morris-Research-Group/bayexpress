#!/usr/bin/env Rscript

# Quick script to do DGE analysis using DEseq2, can be used independently but meant to be called by a python script 

args = commandArgs(trailingOnly=TRUE)

library(DESeq2)


DESeq_DGE <- function(counts_file, metadata_file, output_file_name){
  countData <- read.csv(counts_file, header = TRUE, sep = ",") 
  metaData <- read.csv(metadata_file, header = TRUE, sep = ",")
  
  dds <- DESeqDataSetFromMatrix(countData=countData, 
                                colData=metaData, 
                                design=~dex, tidy = TRUE)
  dds <- DESeq(dds)
  res <- results(dds)
  write.csv(as.data.frame(res), 
            file=output_file_name)
  
}

DESeq_DGE(args[1], args[2], args[3])

















