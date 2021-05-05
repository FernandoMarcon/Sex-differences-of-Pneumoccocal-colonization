rm(list = ls())
pkgs <- c('tidyverse','BiocParallel','DESeq2','pheatmap','RColorBrewer')
suppressPackageStartupMessages(sapply(pkgs, require, character.only = T))
register(MulticoreParam(4))
group.names <- c('POS_M','POS_F','NEG_M','NEG_F')
dataset.names <- c('Adults1','Adults2','Adults3','Elderly1')
gtm.db = 'KEGG_2019'

#### =============== FUNCTIONS ================ ####
plotPathway <- function(selected.pathway) {
  nes.l <- nes %>% rownames_to_column('pathways') %>% gather('volunteer_id', 'nes',-pathways) %>% filter(pathways == selected.pathway)
  nes.l <- pheno %>% rownames_to_column('volunteer_id') %>%
    merge(.,nes.l, by = 'volunteer_id', all.y = T)
  ggplot(nes.l, aes(dataset, nes, col = dataset)) + geom_boxplot(show.legend = F) + geom_jitter(show.legend = F, size = 3) +
    facet_grid(.~group) + theme_linedraw() + labs(title = selected.pathway)
}

#### =============== ANALYSE PATHWAYS =============== ####
nes_padj <- read.delim('intermediate/volunteer_wise_analysis/ssGSEA/NES_padj0.1_logFC_ssGSEAinput.csv', row.names = 2)[,-1]
colnames(nes_padj) = gsub('\\.','\\/',gsub('X','',colnames(nes_padj)))
nes_padj = cbind(pathway = rownames(nes_padj),as.data.frame(apply(nes_padj, 2, as.numeric))) %>%
  column_to_rownames('pathway')
nes_padj[is.na(nes_padj)] <- 0

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

pdf(paste0('intermediate/volunteer_wise_analysis/ssGSEA/KEGG_2019/NES_selectedPathways_padj_',p.thrs,'_volPerc_',vol.perc*100,'_heatmap.pdf'))
plt.nes
dev.off()

#--- NES Boxplot by group
head(pheno)
selected.pathway = selected.pathways[1]
plotPathway(selected.pathway)

pdf(paste0('intermediate/volunteer_wise_analysis/ssGSEA/KEGG_2019/NES_selectedPathways_padj_',p.thrs,'_volPerc_',vol.perc*100,'_boxplot.pdf'))
lapply(selected.pathways, plotPathway)
dev.off()

#---
nes <- nes[selected.pathways,]
nes.l <- nes %>% rownames_to_column('pathway') %>% gather('volunteer_id','nes', -pathway)
nes.l <- pheno %>% rownames_to_column('volunteer_id') %>% unite('class',dataset, group) %>% select(-carriage, -sex) %>%
  merge(.,nes.l, by = 'volunteer_id')

# nes.mean <- nes.l %>% group_by(class, pathway) %>% summarize(nes_mean = median(nes))
# temp = nes.mean %>% spread(class, nes_mean) %>% column_to_rownames('pathway')
# pheatmap(temp, show_rownames = F, col_annotation = pheno)
nes.l <- nes.l %>% mutate(class = gsub('Adults1|Adults2|Adults3','Adults',nes.l$class),
                          nes_abs = abs(nes))
anova_one_way <- lapply(split(nes.l, nes.l$pathway), aov, formula = nes_abs~class)
anova.pval <- sapply(anova_one_way, function(x) unlist(summary(x))[["Pr(>F)1"]])
anova.pathways <- sort(anova.pval[which(anova.pval < .1)]) %>% as.data.frame %>%
  setNames('padj') %>% rownames_to_column('pathway') %>%
  mutate(padj = -log10(padj))

plt.anovaPval <- ggplot(anova.pathways, aes(padj, reorder(pathway,padj), fill = padj)) +
  geom_bar(stat = 'identity') +
  labs(x = '-log10(padj)', title = 'ANOVA Test',subtitle = gtm.db, y = '')

pdf(file.path('intermediate/volunteer_wise_analysis/ssGSEA',gtm.db,paste0('anova_selected_pathways_',gtm.db,'.pdf')))
plt.anovaPval
dev.off()

temp = nes.l[which(nes.l$pathway %in% anova.pathways$pathway),] %>%
  separate(class, c('dataset', 'carriage','sex'), sep ='_')


pdf(file.path('intermediate/volunteer_wise_analysis/ssGSEA',gtm.db,paste0('anova_selected_pathways_',gtm.db,'_boxplot.pdf')))
lapply(unique(temp$pathway),function(selected.pathway) {
  temp %>% filter(pathway == selected.pathway) %>%
    ggplot(aes(dataset, nes, fill = dataset, color = dataset)) +
    geom_jitter(show.legend = F) + geom_boxplot(show.legend = F, alpha = .3) +
      facet_grid(.~sex+carriage) + theme_linedraw() + labs(x = '', y = 'NES',title= selected.pathway) +
      theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank())

  })
dev.off()
