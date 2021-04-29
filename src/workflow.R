rm(list = ls())
pkgs <- c('tidyverse','BiocParallel','pheatmap','RColorBrewer')
suppressPackageStartupMessages(sapply(pkgs, require, character.only = T))
register(MulticoreParam(4))
group.names <- c('POS_M','POS_F','NEG_M','NEG_F')
dataset.names <- c('Adults1','Adults2','Adults3','Elderly1')

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

  pheno <- pheno %>% separate(class, c('dataset','carriage','sex')) %>% unite('class',dataset, sex)

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
  # keep <- edgeR::filterByExpr(de)
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

calcPcombined <- function(df) {
  Pcombined <- unlist(lapply(rownames(df), function(x) {
    metap::sumlog(as.numeric(df[x,])[!is.na(as.numeric(df[x,]))])$p
    }))
  df$FisherMethodP <- Pcombined
  return(df)
}


#### =============== RUN ============== ####
#: DESeq2
res.all <- lapply(dataset.names, function(dataset.name){ # dataset.name = 'Adults1'
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
              # temp <- Reduce(rbind, res.all) %>% filter(pvalue < 0.01, abs(log2FoldChange) > 0)%>% select(gene_id, log2FoldChange, class) %>% separate(class, c('dataset', 'carriage', 'sex')) %>%
              #   unite('class',carriage, sex, sep = '_') %>% mutate(deg = ifelse(log2FoldChange > 0, 'Up', 'Down')) %>% group_by(gene_id, deg, class) %>% summarize(num_studies = n()) %>%
              #   group_by(class, num_studies, deg) %>% summarize(num_degs = n()) %>% group_by(class, num_studies) %>% mutate(total_num_DEGs = sum(num_degs)) %>%
              #   spread(deg, num_degs, fill = 0)
              # temp
              #
              # temp <- Reduce(rbind, res.all) %>% filter(padj < 0.01, abs(log2FoldChange) > 0)%>% select(gene_id, log2FoldChange, class) %>% separate(class, c('dataset', 'carriage', 'sex')) %>%
              #   unite('class',carriage, sex, sep = '_') %>% mutate(deg = ifelse(log2FoldChange > 0, 'Up', 'Down')) %>% group_by(gene_id, deg, class) %>% summarize(num_studies = n()) %>%
              #   group_by(class, num_studies, deg) %>% summarize(num_degs = n()) %>% group_by(class, num_studies) %>% mutate(total_num_DEGs = sum(num_degs)) %>%
              #   spread(deg, num_degs, fill = 0)
              # temp
#--- logFC correlation
# make table with logFCs
logFC.df <- lapply(names(res.all), function(class.name) res.all[[class.name]][,c('gene_id','log2FoldChange')] %>%setNames(c('gene_id',class.name))) %>%
  Reduce(function(x,y) merge(x,y,by = 'gene_id'),.) %>%
    column_to_rownames('gene_id')
# logFC.df <- logFC.df[,grep('Adults',colnames(logFC.df))]
head(logFC.df)

with(logFC.df,plot(Adults1_F, Adults2_F))
with(logFC.df,plot(Adults1_F, Adults3_F))
with(logFC.df,plot(Adults2_F, Adults3_F))
with(logFC.df,plot(Adults1_M, Adults2_M))
with(logFC.df,plot(Adults1_M, Adults3_M))
with(logFC.df,plot(Adults2_M, Adults3_M))
# filter genes
perc.genes = .5
hm_breaks <- seq(-1, 1, length.out = 100)
hm_color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(hm_breaks))
cor.method = 'spearman'

# maxMean
logFC.mean = sort(abs(rowMeans(logFC.df)), decreasing = T)
selected_genes = names(logFC.mean[1:(nrow(logFC.df)*perc.genes)])
plt.maxMean <- pheatmap::pheatmap(cor(logFC.df[selected_genes,], method = cor.method), legend_labels = 'cor',
  main = paste0('maxMean (',perc.genes*100,'% of top genes)'),  color = hm_color,  breaks = hm_breaks)

pdf(file.path('intermediate/item6_DEtables/',paste0('corrplot_logFC_all_maxMean_',gsub('.*\\.','p',as.character(perc.genes)),'.pdf')))
plt.maxMean
dev.off()

# maxVar
logFC.var = sort(apply(logFC.df, 1, var), decreasing = T)
selected_genes = names(logFC.var[1:(nrow(logFC.df)*.1)])
plt.maxVar <- pheatmap::pheatmap(cor(logFC.df[selected_genes,], method = cor.method), legend_labels = 'cor',main = paste0('maxVar (',perc.genes*100,'% of top genes)'),  color = hm_color,  breaks = hm_breaks)

pdf(file.path('intermediate/item6_DEtables/',paste0('corrplot_logFC_all_maxVar_',gsub('.*\\.','p',as.character(perc.genes)),'.pdf')))
plt.maxVar
dev.off()

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
runPcombined_byGroup <- function(group.names, pvalues.df){
  lapply(group.names, function(group.name){
    pvalues.df[, grep(group.name, colnames(pvalues.df))] %>%
      calcPcombined %>% select(FisherMethodP) %>% setNames(group.name) %>% rownames_to_column('gene_id')
    }) %>% Reduce(function(x,y) merge(x,y,by = 'gene_id',all = T),.) %>% column_to_rownames('gene_id')
  }
Pcombined <- runPcombined_byGroup(group.names, pvalues.df)
head(Pcombined)
apply(Pcombined < 0.01,2, sum, na.rm=T)
