rm(list = ls())
pkgs <- c('tidyverse','pheatmap','RColorBrewer')
suppressPackageStartupMessages(sapply(pkgs, require, character.only = T))
group.names <- c('POS_M','POS_F','NEG_M','NEG_F')
dataset.names <- c('Adults1','Adults2','Adults3','Elderly1')
basedir = 'intermediate/volunteer_wise_analysis'
outdir <- file.path(basedir, 'EDA')
if(!dir.exists(outdir)) dir.create(outdir)


#### ============================================= EDA ============================================= ####
logFC.df <- read.delim(file.path(basedir, 'logFC_data.csv'), row.names = 'gene_id')
colnames(logFC.df) <- gsub('X','',gsub('\\.','\\/',colnames(logFC.df)))
pheno <- read.delim(file.path(basedir, 'logFC_pheno.csv'))
rownames(pheno) <- pheno$volunteer_id

pcr <- data.frame(prcomp(t(logFC.df))$x)
pcr$class <- pheno$class
pcr = pcr %>% separate(class, c('dataset','carriage','sex'), sep = '_', remove = F)
pcr = pcr %>% unite('group', carriage, sex, sep = '_', remove = F)

pdf(file.path(outdir, 'pca_all.pdf'))
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
annotation_row = pheno %>% separate(class, c('dataset','carriage','sex')) %>% select(-volunteer_id, -dataset) %>% unite('class', carriage, sex, sep = '_')

# maxMean
logFC.mean = sort(abs(rowMeans(logFC.df)), decreasing = T)
selected_genes = names(logFC.mean[1:(nrow(logFC.df)*perc.genes)])
cor.mtx <- cor(logFC.df[selected_genes,], method = cor.method)

plt.maxMean <- pheatmap(cor.mtx, legend_labels = 'cor', annotation_col = annotation_col,annotation_row = annotation_row,
  show_rownames = F, show_colnames = F,  color = hm_color,  breaks = hm_breaks, main = 'Correlation between volunteers logFC (maxMean)')

pdf(file.path(outdir,paste0('vol_logFC_corSpearman_maxMean_',gsub('.*\\.','p',as.character(perc.genes)),'.pdf')))
plt.maxMean
dev.off()

# maxMean
logFC.var = sort(apply(logFC.df, 1, var), decreasing = T)
selected_genes = names(logFC.var[1:(nrow(logFC.df)*.1)])
cor.mtx <- cor(logFC.df[selected_genes,], method = cor.method)

plt.maxVar <- pheatmap(cor.mtx, legend_labels = 'cor', annotation_col = annotation_col,annotation_row = annotation_row,
  show_rownames = F, show_colnames = F,  color = hm_color,  breaks = hm_breaks, main = 'Correlation between volunteers logFC (maxVar)')
pdf(file.path(outdir,paste0('vol_logFC_corSpearman_maxVar_',gsub('.*\\.','p',as.character(perc.genes)),'.pdf')))
plt.maxVar
dev.off()
