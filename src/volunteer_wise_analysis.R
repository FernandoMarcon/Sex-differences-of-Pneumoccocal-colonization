rm(list = ls())
pkgs <- c('tidyverse','BiocParallel','DESeq2','pheatmap','RColorBrewer')
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

  # pheno <- pheno %>% separate(class, c('dataset','carriage','sex')) %>% unite('class',dataset, sex)

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

#### =============== RUN =============== ####
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
identical(rownames(pheno), colnames(logFC.df))

pcr <- data.frame(prcomp(t(logFC.df))$x)
pcr$class <- pheno$class
pcr = pcr %>% separate(class, c('dataset','carriage','sex'), sep = '_', remove = F)
pcr = pcr %>% unite('group', carriage, sex, sep = '_', remove = F)

pdf('intermediate/volunteer_wise_analysis/pca_all.pdf')
ggplot(pcr, aes(PC1, PC2, col= group)) + geom_point(size = 2) + theme_minimal() + theme(legend.position = 'top')
ggplot(pcr, aes(PC1, PC2, col= dataset)) + geom_point(size = 2) + theme_minimal() + theme(legend.position = 'top')
ggplot(pcr, aes(PC1, PC2, col= carriage)) + geom_point(size = 2) + theme_minimal() + theme(legend.position = 'top')
ggplot(pcr, aes(PC1, PC2, col= sex)) + geom_point(size = 2) + theme_minimal() + theme(legend.position = 'top')
ggplot(pcr, aes(PC1, PC2, col= class)) + geom_point(size = 2) + theme_minimal() + theme(legend.position = 'top')
dev.off()
# filter genes
perc.genes = .5
hm_breaks <- seq(-1, 1, length.out = 100)
hm_color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(hm_breaks))
cor.method = 'spearman'
annotation_col = pheno %>% separate(class, c('dataset','carriage','sex')) %>% select(-volunteer_id)
annotation_row = pheno %>% mutate(class = gsub('Adults|Elderly','',class))%>% select(-volunteer_id)

# maxMean
logFC.mean = sort(abs(rowMeans(logFC.df)), decreasing = T)
selected_genes = names(logFC.mean[1:(nrow(logFC.df)*perc.genes)])
cor.mtx <- cor(logFC.df[selected_genes,], method = cor.method)

plt.maxMean <- pheatmap(cor.mtx, legend_labels = 'cor', annotation_col = annotation_col,annotation_row = annotation_row,
  show_rownames = F, show_colnames = F,  color = hm_color,  breaks = hm_breaks, main = 'Correlation between volunteers logFC (maxMean)')
pdf(file.path('intermediate/volunteer_wise_analysis/',paste0('vol_logFC_corSpearman_maxMean_',gsub('.*\\.','p',as.character(perc.genes)),'.pdf')))
plt.maxMean
dev.off()

# maxMean
logFC.var = sort(apply(logFC.df, 1, var), decreasing = T)
selected_genes = names(logFC.var[1:(nrow(logFC.df)*.1)])
cor.mtx <- cor(logFC.df[selected_genes,], method = cor.method)

plt.maxVar <- pheatmap(cor.mtx, legend_labels = 'cor', annotation_col = annotation_col,annotation_row = annotation_row,
  show_rownames = F, show_colnames = F,  color = hm_color,  breaks = hm_breaks, main = 'Correlation between volunteers logFC (maxVar)')
pdf(file.path('intermediate/volunteer_wise_analysis/',paste0('vol_logFC_corSpearman_maxVar_',gsub('.*\\.','p',as.character(perc.genes)),'.pdf')))
plt.maxVar
dev.off()

# Save data
logFC.df %>% rownames_to_column('gene_id') %>%
  write.table('intermediate/volunteer_wise_analysis/logFC.csv', sep = '\t',row.names = F, quote = F)
write.table(pheno, 'intermediate/volunteer_wise_analysis/logFC_pheno.csv', sep = '\t',row.names = F, quote = F)

#### =============== ssGSEA =============== ####
dir.create('intermediate/volunteer_wise_analysis/ssGSEA', showWarnings = F)

logFC.df = read.delim('intermediate/volunteer_wise_analysis/logFC.csv')
colnames(logFC.df) = gsub('X','',gsub('\\.','\\/',colnames(logFC.df)))

pheno <- read.delim('intermediate/volunteer_wise_analysis/logFC_pheno.csv')
rownames(pheno) <- pheno$volunteer_id

# ensembl to gene symbol
gene.dic <- read.delim('data/tidy_data/gene_annotation.csv')
temp = merge(gene.dic, logFC.df, by = 'gene_id')

# remove duplicated genes
temp = temp %>% filter(symbol != '') %>%
  mutate(gene_mean = rowMeans(.[,-c(1,2)])) %>%
  group_by(symbol) %>% filter(gene_mean == max(gene_mean)) %>%
  select(-gene_mean) %>% select(-gene_id) %>%
  rename(symbol = 'genes')
head(temp)
write.table(temp, 'intermediate/volunteer_wise_analysis/ssGSEA/logFC_ssGSEAinput.csv', sep = '\t',row.names = F, quote = F)

basedir <- '/home/marcon/Documents/work/Sex-differences-of-Pneumoccocal-colonization/'
setwd(file.path(basedir,'intermediate/volunteer_wise_analysis/ssGSEA/'))

# run Single_Sample_GSEA_ssGSEA_fgsea.R
source(file.path(basedir,'src/Single_Sample_GSEA_ssGSEA_fgsea.R'))
gmtfile <- file.path(basedir,'data/Reactome_2016')
fileranks <- "logFC_ssGSEAinput.csv"
Ptype <- "padj"
pval_cutoff <- 0.1
ssGSEA(gmtfile=gmtfile,fileranks=fileranks,Ptype=Ptype,pval_cutoff=pval_cutoff)

#### =============== ANALYSE PATHWAYS =============== ####
# nes_padj <- read.delim('NES_logFC_ssGSEAinput.csv', row.names = 1)
nes_padj <- read.delim('intermediate/volunteer_wise_analysis/ssGSEA/NES_padj0.1_logFC_ssGSEAinput.csv', row.names = 2)[,-1]
colnames(nes_padj) = gsub('\\.','\\/',gsub('X','',colnames(nes_padj)))
nes_padj = cbind(pathway = rownames(nes_padj),as.data.frame(apply(nes_padj, 2, as.numeric))) %>%
  column_to_rownames('pathway')
nes_padj[is.na(nes_padj)] <- 0
# nes_padj = nes_padj %>% mutate(row_mean = rowMeans(.)) %>%
  # filter(abs(row_mean) > 0)

pheno <- read.delim('intermediate/volunteer_wise_analysis/logFC_pheno.csv', row.names = 1) %>%
  separate('class',c('dataset','carriage', 'sex')) %>%
  unite('group', carriage, sex, remove = F)

plt.nes <- pheatmap(nes_padj, show_rownames = F, show_colnames = F, annotation_col = pheno)

pdf('intermediate/volunteer_wise_analysis/ssGSEA/NES_padj0.1_logFC_allPaths_heatmap.pdf')
plt.nes
dev.off()
#--- Filter Pathways
# Filter by padj < threshold for at least X% of volunteers in a given group
nes <- read.delim('intermediate/volunteer_wise_analysis/ssGSEA/NES_logFC_ssGSEAinput.csv', row.names = 1)
padj <- read.delim('intermediate/volunteer_wise_analysis/ssGSEA/padj_logFC_ssGSEAinput.csv', row.names = 1)
colnames(padj) <- colnames(nes) <- gsub('\\.','\\/',gsub('X','',colnames(nes)))
pathways = rownames(nes)
nes <- as.data.frame(apply(nes, 2, as.numeric))
padj <- as.data.frame(apply(padj, 2, as.numeric))
rownames(nes) <- rownames(padj) <- pathways
identical(dimnames(nes),dimnames(padj))

p.thrs = 0.01
vol.perc = .7
selected.pathways <- lapply(group.names, function(group.name){ #group.name = group.names[1]
  selected.vol = pheno %>% filter(group == group.name) %>% rownames
  temp <- padj[,intersect(colnames(padj), selected.vol)]
  selected.pathways <- names(rowSums(temp < p.thrs) > ncol(temp)*vol.perc)
  selected.pathways
  }) %>% unlist %>% unique

plt.nes <- pheatmap(nes[selected.pathways, ], show_rownames = F, show_colnames = F, annotation_col = pheno,
  main = paste0('Pathways with padj < ',p.thrs,' for > ',vol.perc*100,'% of volunteers in a group'))

pdf(paste0('intermediate/volunteer_wise_analysis/ssGSEA/NES_selectedPathways_padj_',p.thrs,'_volPerc_',vol.perc*100,'_heatmap.pdf'))
plt.nes
dev.off()

#--- NES Boxplot by group
boxplot(t(nes[selected.pathways,]))


head(nes.l)

plotPathway <- function(selected.pathway) {
  nes.l <- nes %>% rownames_to_column('pathways') %>% gather('volunteer_id', 'nes',-pathways) %>% filter(pathways == selected.pathway)
  nes.l <- pheno %>% rownames_to_column('volunteer_id') %>%
    merge(.,nes.l, by = 'volunteer_id', all.y = T)
  ggplot(nes.l, aes(sex, nes, col = group)) + geom_boxplot(show.legend = F) + geom_jitter(show.legend = F, size = 3) +
    facet_grid(.~carriage) + theme_linedraw() + labs(title = selected.pathway)
}

selected.pathway = selected.pathways[3]
plotPathway(selected.pathway)

hist(as.numeric(nes[1,]))

nes.l <- nes %>% rownames_to_column('pathways') %>% gather('volunteer_id','nes', -pathways)
nes.l <- pheno %>% rownames_to_column('volunteer_id') %>% merge(., nes.l, by = 'volunteer_id', all.y = T)
head(nes.l)

temp = nes.l %>% filter(pathways == selected.pathway)
res_aov <- aov(nes ~ group, data = temp)

library(car)
par(mfrow = c(1,2))
hist(res_aov$residuals)
qqPlot(res_aov$residuals, id = F)
shapiro.test(res_aov$residuals)
boxplot(nes~group, data = temp)
leveneTest(nes~group, data = temp)


par(mfrow = c(1, 2)) # combine plots
# 1. Homogeneity of variances
plot(res_aov, which = 1)
# 2. Normality
plot(res_aov, which = 2)


#--- Correlation between NES and Spn density
meta <- read.delim('data/raw_data/EHPC_Density_Data_12MARCH2020.tsv', row.names = 1)
meta = meta %>% filter(volunteer_id %in% colnames(nes)) %>%
  mutate(density = ifelse(is.na(density), 0, density), density = log10(density + 1)) %>%
  group_by(volunteer_id) %>% summarize(max_density = max(density))

selected.pathway = selected.pathways[1]
nes.l <- nes %>% rownames_to_column('pathways') %>% gather('volunteer_id','nes', -pathways)
nes.l <- pheno %>% rownames_to_column('volunteer_id') %>% merge(., nes.l, by = 'volunteer_id', all.y = T)
nes.l <- merge(nes.l, meta, by = 'volunteer_id', all.x = T)
head(nes.l)
temp = nes.l %>% filter(pathways == selected.pathway)

ggplot(temp, aes(nes, max_density, col = group)) + geom_point() + theme(legend.position = 'top')
res_lm <- lm(max_density ~ group + nes, data = temp)
names(res_lm)
res_lm$effects

nes.l <- nes %>% rownames_to_column('pathways') %>% gather('volunteer_id','nes', -pathways)
nes.l <- pheno %>% rownames_to_column('volunteer_id') %>% merge(., nes.l, by = 'volunteer_id', all.y = T)
nes.l <- merge(nes.l, meta, by = 'volunteer_id', all.x = T)

cor.res <- lapply(selected.pathways, function(selected.pathway){
  temp = nes.l %>% filter(pathways == selected.pathway)
  return(with(temp, cor(max_density, nes, method = 'spearman')))
  }) %>% unlist %>% setNames(selected.pathways)

boxplot(cor.res)


density.paths <- names(which(abs(cor.res) > .2))
plt.nes <- pheatmap(nes[density.paths, ], show_rownames = F, show_colnames = F, annotation_col = pheno,
  main = paste0('Pathways with padj < ',p.thrs,' for > ',vol.perc*100,'% of volunteers in a group'))

#---- Multivariate Analysis on Pahtways
library(Factoshiny)
library(FactoMineR)

meta <- read.delim('data/raw_data/EHPC_Density_Data_12MARCH2020.tsv', row.names = 1)
meta = meta %>% filter(volunteer_id %in% colnames(nes)) %>%
  mutate(density = ifelse(is.na(density), 0, density), density = log10(density + 1)) %>%
  group_by(volunteer_id) %>% mutate(max_density = max(density))

mydata <- nes[selected.pathways,] %>% t %>% as.data.frame %>% rownames_to_column('volunteer_id') %>%
  merge(meta,., by = 'volunteer_id',all.y = T) %>% mutate(sample_id = paste0(volunteer_id, '_', timepoint)) %>%
  column_to_rownames('sample_id') %>% select(-study, -vaccine_date, -virus,-timepoint_naturalcarriage,-date_inoculated,
    -group, -age, -vaccine, -virus_species, -serotype_naturalcarriage,-inoculation_dose, -inoculation_serotype
    -volunteer_id) %>%
  mutate(timepoint = as.numeric(gsub('D','',timepoint))) %>%
head(mydata)

# PCAshiny(Mydata)

res.pca = PCA(mydata, quanti.sup=4,quali.sup=c(1:3,5:9))
resshiny = PCAshiny(res.pca)
