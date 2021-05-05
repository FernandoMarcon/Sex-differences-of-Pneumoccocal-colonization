rm(list = ls())
pkgs <- c('tidyverse','BiocParallel','DESeq2')
suppressPackageStartupMessages(sapply(pkgs, require, character.only = T))
register(MulticoreParam(4))
group.names <- c('POS_M','POS_F','NEG_M','NEG_F')
dataset.names <- c('Adults1','Adults2','Adults3','Elderly1')
outdir = 'intermediate/volunteer_wise_analysis'

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

  testConsistency(dataset)
}

generateDesignMtx <- function(dataset) {
  pheno <- dataset$pheno
  dataset$design.mtx <- lapply(levels(pheno$class), function(class.name) {
    data.frame(ifelse(pheno$class == class.name & pheno$timepoint != 'baseline', 1, 0))
    }) %>% Reduce(cbind, .) %>% setNames(levels(pheno$class)) %>% as.matrix %>%
    cbind(model.matrix(~ volunteer_id, data = pheno))
  return(dataset)
  }

normData <- function(dataset) {
  dds <- DESeq2::DESeqDataSetFromMatrix(dataset$counts, dataset$pheno, design = dataset$design)
  dataset$norm_data <- assay(vst(dds))
  return(dataset)
}

calc_log2FC <- function(dataset) {
  pheno <- dataset$pheno %>% rownames_to_column('sample_id')
  lapply(unique(pheno$volunteer), function(volunteer.id) {
    sample.tp <- pheno %>% filter(volunteer_id == volunteer.id, timepoint != 'baseline') %>% .$sample_id
    sample.bl <- pheno %>% filter(volunteer_id == volunteer.id, timepoint == 'baseline') %>% .$sample_id
    data.frame(tp = dataset$norm_data[,sample.tp], bl = dataset$norm_data[,sample.bl]) %>%
      mutate(logFC = log2(tp) - log2(bl)) %>% select(logFC) %>% setNames(volunteer.id) %>%
      rownames_to_column('gene_id')
    }) %>% Reduce(function(x,y) merge(x,y,by = 'gene_id',all = T),.) #%>% column_to_rownames('gene_id')
}

#### ============================================= RUN ============================================= ####
# Create directory
if(!dir.exists(outdir)) dir.create(outdir)

# Generate volunteer-wise logFCs table
logFC.df <- lapply(dataset.names, function(dataset.name) {
  dataset <- loadDataset(dataset.name)
  dataset <- removeUnpairedSamples(dataset)
  dataset <- generateDesignMtx(dataset)
  dataset <- normData(dataset)
  calc_log2FC(dataset)
  }) %>% Reduce(function(x,y) merge(x,y,by = 'gene_id'),.) %>%
    column_to_rownames('gene_id')
head(logFC.df)

pheno <- lapply(dataset.names, function(dataset.name) {
  dataset <- loadDataset(dataset.name)
  dataset$pheno %>% select(volunteer_id, class) %>% filter(!duplicated(volunteer_id))
  }) %>% Reduce(rbind, .)
rownames(pheno) <- pheno$volunteer_id
pheno <- pheno[colnames(logFC.df),]
head(pheno)

# Write out tables
if(identical(rownames(pheno), colnames(logFC.df))){
  logFC.df %>% rownames_to_column('gene_id') %>% write.table(file.path(outdir, 'logFC_data.csv'), sep = '\t', quote = F, row.names = F)
  pheno %>% write.table(file.path(outdir, 'logFC_pheno.csv'), sep = '\t', quote = F, row.names = F)
}
