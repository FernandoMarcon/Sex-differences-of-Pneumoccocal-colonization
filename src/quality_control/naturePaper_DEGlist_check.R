#--- SETTINGS
rm(list = ls())
pkgs <- c('tidyverse', 'GEOquery','DESeq2','BiocParallel','xlsx')
suppressPackageStartupMessages(sapply(pkgs, require, character.only = T))
register(MulticoreParam(4))

#--- DATA
# pheno
pheno <- read.delim('./data/tidy_data/tables/Adults2_pheno.csv', row.names = 1)

# counts
counts <- read.delim('./data/tidy_data/tables/Adults2_counts.csv', row.names = 1)
identical(rownames(pheno), colnames(counts))

#--- RUN
pheno = pheno %>% separate(class, c('dataset','carriage', 'sex'))

# Carriage: POS
pheno.pos <- pheno[which(pheno$carriage == 'POS'),]
counts.pos <- counts[,rownames(pheno.pos)]
identical(rownames(pheno.pos), colnames(counts.pos))
dds.pos <- DESeqDataSetFromMatrix(counts.pos, pheno.pos, design = ~ volunteer_id + timepoint)
dds.pos <- DESeq(dds.pos, parallel = T)
res.pos <- results(dds.pos, name = 'timepoint_D9_vs_baseline', parallel = T)
res.pos <- res.pos %>% as.data.frame %>% select(log2FoldChange, pvalue) %>% rownames_to_column('gene_id')

# Carriage: NEG
pheno.neg <- pheno[which(pheno$carriage == 'NEG'),]
counts.neg <- counts[,rownames(pheno.neg)]
identical(rownames(pheno.neg), colnames(counts.neg))
dds.neg <- DESeqDataSetFromMatrix(counts.neg, pheno.neg, design = ~ volunteer_id + timepoint)
dds.neg <- DESeq(dds.neg, parallel = T)
res.neg <- results(dds.neg, name = 'timepoint_D9_vs_baseline', parallel = T)
res.neg <- res.neg %>% as.data.frame %>% select(log2FoldChange, pvalue) %>% rownames_to_column('gene_id')

# DEG list
# https://static-content.springer.com/esm/art%3A10.1038%2Fs41590-018-0231-y/MediaObjects/41590_2018_231_MOESM6_ESM.xlsx
deg_list <- read.xlsx('data/nature_paper_data/41590_2018_231_MOESM6_ESM.xlsx', 1) %>%
    mutate(Ensembl_ID = gsub('\\..*','',Ensembl_ID)) %>% rename(Ensembl_ID = 'gene_id')

merged <- merge(res.pos[,c('gene_id','log2FoldChange')],
                res.neg[,c('gene_id','log2FoldChange')],
                by = 'gene_id', all = T)
colnames(merged) <- c('gene_id','POS','NEG')
merged <- merge(merged,
                deg_list[,c('gene_id','Control_Carriage._Day9_log2FC','Control_Carriage._Day9_log2FC.1')],
                by = 'gene_id',all = T)
ggplot(merged, aes(POS, Control_Carriage._Day9_log2FC)) + geom_point()
ggplot(merged, aes(NEG, Control_Carriage._Day9_log2FC.1)) + geom_point()
