## Evaluate p-value thresholds

### Settings
rm(list = ls())

basedir <- '~/Documents/work/Sex-differences-of-Pneumoccocal-colonization/'
outdir = file.path(basedir,'intermediate/item6_DEtables/')
if(!dir.exists(outdir)) dir.create(outdir)

pkgs <- c('tidyverse','data.table','ggplot2')
suppressPackageStartupMessages(sapply(pkgs, require, character.only = T))
### Load Data
fnames <- list.files(outdir, full.names = T, pattern = 'DESeq2_topTable')
de.table <- lapply(fnames, fread, data.table = F) %>% Reduce(rbind, .)
head(de.table)

pvalues <- seq(0,.1,.01)

### Non-adjusted p-value
deg_byPvalThrs <- lapply(pvalues, function(pval) {
  de.table %>% filter(abs(log2FoldChange) > 0, pvalue < pval) %>%
    group_by(group) %>%
    summarise(num_deg = n()) %>%
    mutate(pval_thrs = pval)
}) %>% Reduce(rbind, .) %>% separate(group, into = c('dataset','carriage','sex'),remove = F)
head(deg_byPvalThrs)

plt_study <- ggplot(deg_byPvalThrs, aes(-log10(pval_thrs), num_deg, group = group, col = dataset)) +
  geom_point() + geom_line() +  labs(x = '-log10(p-value)', y = 'Number of DEGs') + theme_minimal() +
  geom_hline(yintercept = c(100, 1000), col = 'red') +
  geom_vline(xintercept = c(-log10(0.01),-log10(0.05),-log10(0.1)), linetype ='dashed')
plt_study
# Color by sex
plt_sex <- ggplot(deg_byPvalThrs, aes(-log10(pval_thrs), num_deg, group = group, col = sex)) +
  geom_point() + geom_line() +  labs(x = '-log10(p-value)', y = 'Number of DEGs', col = 'Sex') + theme_minimal() +
  geom_hline(yintercept = c(100, 1000), col = 'red') +
  geom_vline(xintercept = c(-log10(0.01),-log10(0.05),-log10(0.1)), linetype ='dashed')
plt_sex
# Color by carriage
plt_carriage <- ggplot(deg_byPvalThrs, aes(-log10(pval_thrs), num_deg, group = group, col = carriage)) +
  geom_point() + geom_line() +  labs(x = '-log10(p-value)', y = 'Number of DEGs', col = 'Carriage') + theme_minimal() +
  geom_hline(yintercept = c(100, 1000), col = 'red') +
  geom_vline(xintercept = c(-log10(0.01),-log10(0.05),-log10(0.1)), linetype ='dashed')
plt_carriage
# Color by Sex and Carriage
plt_sex_carriage <- ggplot(deg_byPvalThrs, aes(-log10(pval_thrs), num_deg, group = paste0(dataset, group), col = paste0(carriage, '_',sex))) +
  geom_point() + geom_line() +  labs(x = '-log10(p-value)', y = 'Number of DEGs', col = 'Carriage + Sex') + theme_minimal() +
  geom_hline(yintercept = c(100, 1000), col = 'red') +
  geom_vline(xintercept = c(-log10(0.01),-log10(0.05),-log10(0.1)), linetype ='dashed')
plt_sex_carriage
pdf(file.path(outdir, 'numDEGs_by_pval_thresholds_assessment_DESeq2.pdf'))
plt_study
plt_sex
plt_carriage
plt_sex_carriage
dev.off()
### Adjusted p-value
deg_byFDRThrs <- lapply(pvalues, function(pval) {
  de.table %>% filter(abs(log2FoldChange) > 0, padj < pval) %>%
    group_by(group) %>%
    summarise(num_deg = n()) %>%
    mutate(pval_thrs = pval)
}) %>% Reduce(rbind, .) %>% separate(group, into = c('dataset','carriage','sex'),remove = F)

plt.dataset <- ggplot(deg_byFDRThrs, aes(-log10(pval_thrs), num_deg, group = group, col = dataset)) +
  geom_point() + geom_line() + labs(x = '-log10(FDR)', y = 'Number of DEGs', col = 'dataset') + theme_minimal() +
  geom_hline(yintercept = c(100, 1000), col = 'red') +
  geom_vline(xintercept = c(-log10(0.01),-log10(0.05),-log10(0.1)), linetype ='dashed')
plt.dataset

plt.sex <- ggplot(deg_byFDRThrs, aes(-log10(pval_thrs), num_deg, group = group, col = sex)) +
  geom_point() + geom_line() + labs(x = '-log10(FDR)', y = 'Number of DEGs', col = 'sex') + theme_minimal() +
  geom_hline(yintercept = c(100, 1000), col = 'red') +
  geom_vline(xintercept = c(-log10(0.01),-log10(0.05),-log10(0.1)), linetype ='dashed')
plt.sex

plt.carriage <- ggplot(deg_byFDRThrs, aes(-log10(pval_thrs), num_deg, group = group, col = carriage)) +
  geom_point() + geom_line() + labs(x = '-log10(FDR)', y = 'Number of DEGs', col = 'carriage') + theme_minimal() +
  geom_hline(yintercept = c(100, 1000), col = 'red') +
  geom_vline(xintercept = c(-log10(0.01),-log10(0.05),-log10(0.1)), linetype ='dashed')
plt.carriage

plt.group <- ggplot(deg_byFDRThrs, aes(-log10(pval_thrs), num_deg, group = group, col = paste0(carriage, '_', sex))) +
  geom_point() + geom_line() + labs(x = '-log10(FDR)', y = 'Number of DEGs', col = 'group') + theme_minimal() +
  geom_hline(yintercept = c(100, 1000), col = 'red') +
  geom_vline(xintercept = c(-log10(0.01),-log10(0.05),-log10(0.1)), linetype ='dashed')
plt.group

pdf(file.path(outdir, 'numDEGs_by_FDR_thresholds_assessment_DESeq2.pdf'))
plt.dataset
plt.sex
plt.carriage
plt.group
dev.off()

## Plot Number of DEGs by Group
pval = 0.01
num_degs <- de.table %>% filter(pvalue < pval) %>%
  mutate(deg = ifelse(log2FoldChange > 0,1, 0), deg = ifelse(log2FoldChange < 0, -1, deg)) %>%
  group_by(group, deg) %>%
  summarise(num_deg = n()) %>%
  separate(group, into = c('dataset','carriage','sex'),remove = F)

plt.num_degs <- ggplot(num_degs, aes(sex, deg*num_deg, fill = ifelse(deg > 0, 'Up', 'Down'))) +
  geom_bar(stat = 'identity') +
  facet_grid(.~dataset + carriage) +theme_linedraw() +
  scale_fill_manual(values = c('dodgerblue3','darkred')) +
  labs(x = 'Sex', y = 'Number of DEGs', fill = 'Direction') +
  theme(panel.grid = element_blank())
plt.num_degs

pdf(file.path(outdir, 'num_degs_pval-01_barplot_DESeq2.pdf'), width = 10)
plt.num_degs
dev.off()

# Save Num. DEGs Table
num_degs %>%
  mutate(deg = ifelse(deg == -1, 'Down', 'Up')) %>%
  spread(key = deg,value = num_deg) %>%
  fwrite(file.path(outdir, 'num_deg_pval-01_table.csv'), sep = '\t', quote = F, row.names = F)
