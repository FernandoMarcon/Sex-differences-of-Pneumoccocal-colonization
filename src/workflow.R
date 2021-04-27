rm(list = ls())
pkgs <- c('tidyverse','BiocParallel')
suppressPackageStartupMessages(sapply(pkgs, require, character.only = T))
register(MulticoreParam(4))

#### =============== FUNCTIONS =============== ####
testConsistency <- function(dataset) {
  if(identical(rownames(dataset$pheno), colnames(dataset$counts))) {
    return(dataset)
  } else {
    stop('Names not identical')
  }
}

loadDataset <- function(dataset.name) {
  pheno <- read.delim(file.path('./data/tidy_data/tables/',paste0(dataset.name, '_pheno.csv')), row.names = 'sample_id')
  counts <- read.delim(file.path('./data/tidy_data/tables/',paste0(dataset.name, '_counts.csv')), row.names = 'gene_id')
  testConsistency(list(pheno = pheno, counts = counts))
}

removeUnpairedSamples <- function(dataset) {
  pheno <- dataset$pheno
  num_samples = pheno %>% group_by(volunteer_id) %>% summarise(num_samples = n())
  selected_vol = num_samples %>% filter(num_samples > 1) %>% .$volunteer_id
  dataset$pheno <- pheno[which(pheno$volunteer_id %in% selected_vol),]
  dataset$counts <- dataset$counts[,rownames(dataset$pheno)]
  message('Removed:\n')
  num_samples %>% filter(num_samples == 1) %>% print
  testConsistency(dataset)
}

filterCounts <- function(dataset) {
  de <- edgeR::DGEList(counts = dataset$counts, samples = dataset$pheno)
  keep <- edgeR::filterByExpr(de, group = 'class')
  dataset$counts <- dataset$counts[keep, ]
  message("\nRemoved ",sum(keep),' genes, remaining: ',sum(!keep),'\n')
  return(dataset)
}

dataset.name = 'Adults1'
dataset <- loadDataset(dataset.name)
dataset <- removeUnpairedSamples(dataset)
dataset <- filterCounts(dataset)
runDE <- function(dataset, de.method = c('edgeR','DESeq2')){

}

runDESeq2 <- function(dataset) {
  dds <- DESeq2::DESeqFromMatrix(dataset$counts, dataset$pheno, design = '~ 1')
  return(dds)
}

dds <- runDESeq2(dataset)
