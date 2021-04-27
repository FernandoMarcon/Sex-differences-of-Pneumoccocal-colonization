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

  dataset$pheno$timepoint <- as.factor(dataset$pheno$timepoint) %>% relevel(ref = 'baseline')
  dataset$pheno$class <- as.factor(dataset$pheno$class)
  dataset$pheno$volunteer_id <- as.factor(dataset$pheno$volunteer_id)

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

generateDesignMtx <- function(dataset) {
  pheno <- dataset$pheno
  dataset$design.mtx <- lapply(levels(pheno$class), function(class.name) {
    data.frame(ifelse(pheno$class == class.name & pheno$timepoint != 'baseline', 1, 0))
    }) %>% Reduce(cbind, .) %>% setNames(levels(pheno$class)) %>% as.matrix %>%
    cbind(model.matrix(~ volunteer_id, data = pheno))
  return(dataset)
  # temp = design.mtx %>% as.data.frame %>% filter(Adults1_NEG_F == 1) %>% row.names
  # table(pheno[temp,c('class', 'timepoint')])
  }

runDE <- function(dataset, FUN = ...){
  res <- FUN(dataset)
  return(res)
}

runDESeq2 <- function(dataset) {
  dds <- DESeq2::DESeqDataSetFromMatrix(dataset$counts, dataset$pheno, design = dataset$design)
  dds <- DESeq2::DESeq(dds, parallel = T)
  de.res <- lapply(levels(dataset$pheno$class), function(class.name) {
    res <- results(dds, name = class.name)
    res %>% as.data.frame %>% mutate(class = class.name) %>%
    rownames_to_column('gene_id') %>%
    rename(log2FoldChange = 'logFC')
  }) %>% setNames(levels(dataset$pheno$class))
  return(de.res)
}

#--- RUN: DESeq2
datasets <- c('Adults1','Adults2','Adults3','Elderly1')
res.all <- lapply(datasets, function(dataset.name){
  # dataset.name = 'Adults1'
  message('\n\n',dataset.name, '\n\n')
  dataset <- loadDataset(dataset.name)
  dataset <- removeUnpairedSamples(dataset)
  dataset <- filterCounts(dataset)
  dataset <- generateDesignMtx(dataset)
  res <- runDE(dataset, runDESeq2)
  return(res)
  }) %>% unlist(recursive = F)


names(res.all)

#--- Pcombined (Fischer method)
#--- MetaVolcano (vote count)
