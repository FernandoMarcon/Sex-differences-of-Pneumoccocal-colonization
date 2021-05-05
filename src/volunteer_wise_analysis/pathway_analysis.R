rm(list = ls())
pkgs <- c('tidyverse','pheatmap','RColorBrewer')
suppressPackageStartupMessages(sapply(pkgs, require, character.only = T))
group.names <- c('POS_M','POS_F','NEG_M','NEG_F')
# dataset.names <- c('Adults1','Adults2','Adults3','Elderly1')

gtm.dbs <- list.files('data/GTM')
#### =============== FUNCTIONS ================ ####
plotPathway <- function(selected.pathway) {
  nes.l <- nes %>% rownames_to_column('pathways') %>% gather('volunteer_id', 'nes',-pathways) %>% filter(pathways == selected.pathway)
  nes.l <- pheno %>% rownames_to_column('volunteer_id') %>%
    merge(.,nes.l, by = 'volunteer_id', all.y = T)
  ggplot(nes.l, aes(dataset, nes, col = dataset)) + geom_boxplot(show.legend = F) + geom_jitter(show.legend = F, size = 3) +
    facet_grid(.~group) + theme_linedraw() + labs(title = selected.pathway)
}

#--- SELECT DB
gtm.db = gtm.dbs[1]
basedir <- file.path('intermediate/volunteer_wise_analysis/ssGSEA',gtm.db)
list.files(basedir)

#--- LOAD DATA
pheno <- read.delim(file.path('intermediate/volunteer_wise_analysis/logFC_pheno.csv'), row.names = 'volunteer_id')
pheno <- pheno %>% separate(class, c('dataset', 'carriage','sex'), sep = '_') %>% unite('group', carriage, sex, sep = '_', remove = F)

nes <- read.delim(file.path(basedir, 'NES_logFC_geneClean.csv'), row.names = 1)
padj <- read.delim(file.path(basedir, 'padj_logFC_geneClean.csv'), row.names = 1)
colnames(padj) <- colnames(nes) <- gsub('X','',gsub('\\.','\\/',colnames(nes)))

#--- FILTERING
# Filter by padj < threshold for at least X% of volunteers in a given group
p.thrs = 0.01
vol.perc = .7
selected.pathways <- lapply(group.names, function(group.name){ #group.name = group.names[1]
  selected.vol = pheno %>% filter(group == group.name) %>% rownames
  temp <- padj[,intersect(colnames(padj), selected.vol)]
  selected.pathways <- names(rowSums(temp < p.thrs) > ncol(temp)*vol.perc)
  selected.pathways
  }) %>% unlist %>% unique

nes = nes[selected.pathways, ]
plt.nes <- pheatmap(nes, show_rownames = F, show_colnames = F, annotation_col = pheno,
  main = paste0('Pathways with padj < ',p.thrs,' for > ',vol.perc*100,'% of volunteers in a group'))

pdf(file.path(basedir, paste0('NES_selectedPathways_padj_',p.thrs,'_volPerc_',vol.perc*100,'_heatmap.pdf')))
plt.nes
dev.off()



head(pheno)
female.vols <- pheno %>% filter(sex == 'F') %>% rownames
male.vols <- pheno %>% filter(sex == 'M') %>% rownames

nes[,female.vols]

tTest_res <- apply(nes, 1, function(x) t.test(abs(x[female.vols]), abs(x[male.vols]))$p.value) %>%
  as.data.frame %>% setNames('pval') %>% rownames_to_column('pathway')
sum(tTest_res$pval < 0.01)
