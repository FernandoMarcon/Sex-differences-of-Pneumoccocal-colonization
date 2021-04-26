rm(list = ls())
pkgs <- c('tidyverse','data.table','ggplot2','DESeq2','BiocParallel','edgeR',"pheatmap","vsn",'RColorBrewer')
suppressPackageStartupMessages(sapply(pkgs, require, character.only = T))
register(MulticoreParam(4))

input.dir <- 'data//tidy_data//tables/'
datasets <- c('Adults1','Adults2','Adults3','Elderly1')

dataset.name = datasets[4]
message(dataset.name)
counts <- fread(file.path(input.dir, paste0(dataset.name, '_counts.csv'))) %>% data.frame(row.names = 1)
pheno <- fread(file.path(input.dir, paste0(dataset.name, '_pheno.csv'))) %>% data.frame(row.names = 1)

# Remove incomplete volunteers (only one timepoint)
vol_keep <- pheno %>% group_by(volunteer_id) %>% summarize(num_samples = n()) %>% filter(num_samples > 1) %>% .$volunteer_id
pheno <- pheno[which(pheno$volunteer_id %in% vol_keep),]
counts <- counts[,rownames(pheno)]

# Groups as factors
pheno$volunteer_id <- as.factor(pheno$volunteer_id)
pheno$timepoint <- as.factor(pheno$timepoint) %>% relevel(ref = 'baseline')
pheno$class <- as.factor(gsub(paste0(dataset.name,'_'),'',pheno$class))

design <- lapply(levels(pheno$class), function(group.name) {
    data.frame(ifelse(pheno$class == group.name & pheno$timepoint != 'baseline', 1, 0)) %>%
    setNames(group.name)
}) %>% Reduce(cbind, .) %>% as.matrix

design <- cbind(model.matrix(~volunteer_id, data = pheno), design)
head(design[,levels(pheno$class)])

y <- DGEList(counts = counts, samples = pheno)
keep <- filterByExpr(y, design[,levels(pheno$class)])
dds <- DESeqDataSetFromMatrix(counts[keep,], pheno, design = design)
dds <- DESeq(dds, parallel = T)
# dds <- DESeqDataSetFromMatrix(counts, pheno, design = design)
ntd <- normTransform(dds)
meanSdPlot(assay(ntd))

vsd <- vst(dds, blind=T)
meanSdPlot(assay(vsd))

# Data quality assessment by sample clustering and visualization
## Heatmap of the count matrix
select <- order(rowMeans(counts(dds,normalized=TRUE)), decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("class","timepoint")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,cluster_cols=FALSE, annotation_col=df)

pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df)

## Heatmap of the sample-to-sample distances
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$class, vsd$timepoint, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists,
    clustering_distance_cols=sampleDists, col=colors)

## Principal component plot of the samples
plotPCA(vsd, intgroup=c("class", "timepoint"))

# Differential expression analysis
dds.res <- lapply(levels(dds$class), function(group.name) results(dds, name = group.name)) %>% setNames(levels(dds$class))

getDEGs <- function(res) res %>% as.data.frame %>%
  mutate(DEG = ifelse(log2FoldChange > 0 & pvalue < 0.01, 1, 0),DEG = ifelse(log2FoldChange < 0 & pvalue < 0.01, -1, DEG))

countDEGs <- function(res) {
    de.tbl <- getDEGs(res)
    de.tbl %>% filter(DEG != 0) %>% mutate(DEG = ifelse(DEG == 1, 'Up', 'Down')) %>% group_by(DEG) %>% summarise(num_DEGs = n())
}
lapply(dds.res, countDEGs)

## Log fold change shrinkage for visualization and ranking
# resLFC <- lapply(levels(dds$class), function(group.name) lfcShrink(dds, coef = group.name, type="apeglm",parallel = T)) %>% (levels(dds$class))
# lapply(resLFC, countDEGs)
