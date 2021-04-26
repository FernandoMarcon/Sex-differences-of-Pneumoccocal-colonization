# # Item 6 - DEG per cohort
#
# EdgeR use um FDR cutoff que dê entre 100 a 1000 DEGs.
#
# -   Adult1_Positive_Male_D2v0 (pareado)
# -   Adult2_Positive_Male_D2v0 (pareado)
# -   Adult3_Positive_Male_D2v0 (pareado)
# -   Elderly1_Positive_Male_D2v0 (pareado) ...repetir para TODOS que tem pelo menos 2 amostras.

### Settings
rm(list = ls())

basedir <- '~/Documents/work/Sex-differences-of-Pneumoccocal-colonization/'
outdir = file.path(basedir,'intermediate/item6_DEtables/')
if(!dir.exists(outdir)) dir.create(outdir)

pkgs <- c('tidyverse','data.table','edgeR','DESeq2','BiocParallel')
suppressPackageStartupMessages(sapply(pkgs, require, character.only = T))
register(MulticoreParam(4))

# Functions
load_data <- function(dataset.name){
  counts <- read.delim(file.path('./data/tidy_data/tables', paste0(dataset.name, '_counts.csv')), row.names = 1)
  pheno <- read.delim(file.path('./data/tidy_data/tables',paste0(dataset.name, '_pheno.csv')))

  pheno = pheno %>% group_by(volunteer_id) %>%
    mutate(n = n()) %>% filter(n > 1) %>%
    column_to_rownames('sample_id') %>%
    relocate(class, volunteer_id, timepoint) %>%
    mutate(class = gsub(paste0(dataset.name,'_'),'',class)) %>%
    select(-study, -vaccine, -n)

  counts = counts[,rownames(pheno)]
  identical(rownames(pheno), colnames(counts))

  pheno$volunteer_id <- factor(pheno$volunteer_id)
  pheno$class <- factor(pheno$class, levels = unique(pheno$class))
  pheno$timepoint <- as.factor(pheno$timepoint) %>% relevel(ref = 'baseline')

  return(dataset = list(name = dataset.name, pheno = pheno, counts = counts))
}

filter_counts <- function(dataset) {
  de <- edgeR::DGEList(counts = dataset$counts, samples = dataset$pheno)
  keep <- edgeR::filterByExpr(de)
  dataset$counts <- dataset$counts[keep, ]
  return(dataset)
}

create_design <- function(dataset) {
  pheno <- dataset$pheno

  # Define the experimental factors
  design.mtx <- model.matrix(~ volunteer_id, data = pheno)
  groups <- levels(pheno$class)

  experimental_factors <- lapply(groups, function(group){
    data.frame(ifelse(with(pheno, class == group & timepoint != 'baseline') == T, 1, 0))
  }) %>% Reduce(cbind, .) %>% setNames(groups)

  dataset$design <- as.matrix(cbind(design.mtx, experimental_factors))

  return(dataset)
}

run_DESeq2 <- function(dataset){
  dds <- DESeqDataSetFromMatrix(dataset$counts, dataset$pheno, design = dataset$design)
  dds <- DESeq(dds, parallel = T)
  comparisons <- grep('F|M',resultsNames(dds), value=T)
  de.tables <- lapply(comparisons, function(comp.name) {
    results(dds, name = comp.name, parallel = T) %>% as.data.frame %>%
      rownames_to_column('gene_id') %>% mutate(group = paste0(dataset$name,'_', comp.name))
    }) %>% setNames(comparisons)
  return(de.tables)
}

datasets.name = c('Adults1', 'Adults2', 'Adults3', 'Elderly1')
de.res.all <- lapply(datasets.name[1:2], function(dataset.name) {
  # dataset.name = datasets[2]
  message(dataset.name)

  dataset <- load_data(dataset.name)
  dataset <- filter_counts(dataset)
  dataset <- create_design(dataset)
  de.tables <- run_DESeq2(dataset)
  names(de.tables) <- paste0(dataset.name,'_',names(de.tables))
  return(de.tables)
 }) %>% unlist(recursive = F)

names(de.res.all)
lapply(names(de.res.all), function(group_name) de.res.all[[group_name]] %>%
         fwrite(file.path(outdir,paste0('DESeq2_topTable_',group_name,'.csv')), sep = '\t',quote = F, row.names = F))
