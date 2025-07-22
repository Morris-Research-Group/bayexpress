library(tidyverse)
library(edgeR)
library(data.table)

args = commandArgs(trailingOnly=TRUE)

input_file <- "~/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/PhD/Registration_experiments/Shannons_expression_entropy_plots/Stellar/datasets_for_DGE/14vs18_FG.csv"
do_DGE <- function(input_file, output_file_name){
  
  input_data <- read.csv(input_file, header = TRUE, sep = ",")
  count_data <- input_data %>% dplyr::select(-locus_name)
  
  # creating group array? factor? thingy
  group <- factor(c(1,1,1,2,2,2))
  
  # creating edgeR DGEList object
  d0 <- edgeR::DGEList(
    counts = count_data,
    group = group,
    genes = input_data$locus_name
  )
  
  d0 <- calcNormFactors(d0)
  #cutoff <- 0
  #drop <- which(apply(cpm(d0), 1, max) < cutoff)
  #d <- d0[-drop,]
  mm <- model.matrix(~0 + group)  
  
  y <- voom(d0, mm, plot = T)
  
  fit <- lmFit(y, mm)
  contr <- makeContrasts(group1 - group2, levels = colnames(coef(fit)))
  tmp <- contrasts.fit(fit, contr)
  tmp <- eBayes(tmp)
  
  top.table <- topTable(tmp, sort.by = "P", n = Inf)
  
  top.table$Gene <- rownames(top.table)
  write.table(top.table, file = output_file_name, row.names = F, sep = ",", quote = F)
}
  