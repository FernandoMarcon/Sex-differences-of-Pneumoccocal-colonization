rm(list = ls())
pkgs <- c('tidyverse','BiocParallel')
suppressPackageStartupMessages(sapply(pkgs, require, character.only = T))
register(MulticoreParam(4))
group.names <- c('POS_M','POS_F','NEG_M','NEG_F')

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
  # keep <- edgeR::filterByExpr(de, group = 'class')
  keep <- edgeR::filterByExpr(de)
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
    res <- DESeq2::results(dds, name = class.name)
    res %>% as.data.frame %>% mutate(class = class.name) %>%
    rownames_to_column('gene_id')
  }) %>% setNames(levels(dataset$pheno$class))
  return(de.res)
}
# %>%  rename(log2FoldChange = 'logFC')
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

# temp.name = names(res.all)[3]
# for(temp.name in names(res.all)){
#   temp <- read.delim(file.path('./intermediate/item6_DEtables/',paste0('DESeq2_topTable_',temp.name,'.csv')))
#   temp = merge(temp[,c('gene_id','log2FoldChange')],
#               res.all[[temp.name]][,c('gene_id','log2FoldChange')],
#               by = 'gene_id') %>% setNames(c('gene_id','original','workflow'))
#   plt <- ggplot(temp, aes(original, workflow)) + geom_point()
#   print(plt)
#   Sys.sleep(5)
# }

# Common DEGs between cohorts
temp <- Reduce(rbind, res.all)
temp = temp %>% mutate(deg = ifelse(log2FoldChange > 0, 'Up', 'unchanged'), deg = ifelse(log2FoldChange < 0, 'Down', deg))
temp = temp %>% filter(pvalue < 0.01, deg != 0)
temp = temp %>% select(gene_id, class, deg)
temp = temp %>% separate(class, c('dataset','carriage','sex'))
temp = temp %>% unite('class', carriage, sex, sep = '_')

head(temp)
temp %>% group_by(class, deg, gene_id) %>% summarize(num_deg = n()) %>% ggplot(aes(deg, num_deg, fill = deg)) + geom_bar(stat = 'identity') + facet_grid(.~class)






#### =============== META-DEGS =============== ####
#--- Pcombined (Fischer method)
# Select all pvalue columns
pvalues.df <- lapply(names(res.all), function(class.name) {
  df = res.all[[class.name]]
  df = df[,c('gene_id','pvalue')]
  colnames(df) <- c('gene_id',class.name)
  return(df)
  }) %>% Reduce(function(x, y) merge(x, y, by = 'gene_id', all = T),.) %>% column_to_rownames('gene_id')

# Calculate Pcombined for all datasets
Pcombined <- unlist(lapply(rownames(pvalues.df), function(x) {
  metap::sumlog(as.numeric(pvalues.df[x,])[!is.na(as.numeric(pvalues.df[x,]))])$p
  }))
pvalues.df$FisherMethodP <- Pcombined
pvalues.df$FDR <- p.adjust(Pcombined, method = "fdr", n = length(Pcombined))

# Calculate Pcombined for Adults only
comp_adults = grep('Adults', colnames(pvalues.df), value = T)
Pcombined <- unlist(lapply(rownames(pvalues.df[,comp_adults]), function(x) {
  metap::sumlog(as.numeric(pvalues.df[x,comp_adults])[!is.na(as.numeric(pvalues.df[x,comp_adults]))])$p
  }))
pvalues.df$FisherMethodP_adults <- Pcombined
pvalues.df$FDR_adults <- p.adjust(Pcombined, method = "fdr", n = length(Pcombined))

table(pvalues.df$FisherMethodP < 0.01)
table(pvalues.df$FisherMethodP_adults < 0.01)

table(pvalues.df$FDR < 0.01)
table(pvalues.df$FDR_adults < 0.01)

# Select all logFC columns
logFC.df <- lapply(names(res.all), function(class.name) {
  df = res.all[[class.name]]
  df = df[,c('gene_id','log2FoldChange')]
  colnames(df) <- c('gene_id',class.name)
  return(df)
  }) %>% Reduce(function(x, y) merge(x, y, by = 'gene_id', all = T),.) %>% column_to_rownames('gene_id')

deg.adults <- logFC.df[,grep('Adults',colnames(logFC.df))]
deg.adults$FDR_adults <- pvalues.df$FDR_adults
deg.adults = deg.adults %>% filter(FDR_adults < 0.01) %>% select(-FDR_adults)

deg.adults.counts <- lapply(group.names, function(group.name) {
  data.frame(up = rowSums(deg.adults[,grep(group.name, colnames(deg.adults))] > 0),
              down = rowSums(deg.adults[,grep(group.name, colnames(deg.adults))] < 0)) %>%
              setNames(paste0(group.name, '|',colnames(.))) %>%
              rownames_to_column('gene_id')
              }) %>% Reduce(function(x, y) merge(x, y, by = 'gene_id', all = T),.)

deg.adults.counts = deg.adults.counts %>% gather(class, num_degs, -gene_id, na.rm = T)

pdf(file.path('intermediate/item7_metaDEGs','voteCount_Pcombined_Fischer_DESeq2.pdf'))
deg.adults.counts %>% group_by(num_degs,class) %>% summarize(vote_count = n()) %>%
separate(class, c('class','direction'), sep = '\\|') %>%
filter(num_degs > 0) %>%
  ggplot(aes(num_degs, vote_count, fill = class)) + geom_bar(stat = 'identity', show.legend = F) +
  facet_grid(.~class+direction) + theme_minimal() +
  labs(title = 'Vote count', subtitle = 'FDR (Pcombined) < 0.01, |logFC| > 0, DESeq2')
dev.off()


# %>%  separate(class, c('class','direction'), sep = '\\|')
