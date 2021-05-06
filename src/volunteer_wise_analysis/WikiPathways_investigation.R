rm(list = ls())
pkgs <- c('tidyverse')
suppressPackageStartupMessages(sapply(pkgs, require, character.only = T))
gtm.db = 'WikiPathways_2019_Human'
basedir <- file.path('intermediate/volunteer_wise_analysis/ssGSEA',gtm.db)

# LOAD DATA
pheno <- read.delim('~/Documents/work/Sex-differences-of-Pneumoccocal-colonization/intermediate/volunteer_wise_analysis/logFC_pheno.csv')

nes <- read.delim(file.path(basedir, 'NES_logFC_geneClean.csv'), row.names = 1)
colnames(nes) <- gsub('X','',gsub('\\.','\\/',colnames(nes)))

pathways <- read.delim(file.path(basedir, paste0(gtm.db, '_selectedPathways.csv'))) %>%
  filter(pval < 0.01) %>% arrange(pval)

#### =============== FUNCTIONS =============== ####
prepData <- function(pathway) {
  nes[pathway,] %>% t %>% as.data.frame %>% setNames('nes') %>%
    rownames_to_column('volunteer_id') %>%
    merge(pheno, ., by = 'volunteer_id') %>%
    separate(class, c('dataset', 'carriage', 'sex'), sep = '_')
}

boxplotPathway <- function(pathway){
  nes.sub <- prepData(pathway)
  ggplot(nes.sub, aes(sex, nes, col = sex)) +  geom_boxplot(show.legend = F) + geom_jitter(show.legend = F) +
    theme_linedraw() + facet_grid(dataset~carriage) + labs(title = pathway)
}

#--- RUN
boxplotPathway(pathways$pathway[1])
boxplotPathway(pathways$pathway[2])
boxplotPathway(pathways$pathway[3])
boxplotPathway(pathways$pathway[4])
boxplotPathway(pathways$pathway[5])
boxplotPathway(pathways$pathway[6])
boxplotPathway(pathways$pathway[7])
boxplotPathway(pathways$pathway[8])
boxplotPathway(pathways$pathway[9])


le <- read.delim(file.path(basedir, 'LE_logFC_geneClean.csv'))
colnames(le)[1] <- 'pathway'
colnames(le) <- gsub('\\.','\\/',gsub('X', '', colnames(le)))
head(le)
le %>% filter(pathway == pathways$pathway[1]) %>% .[,-1] %>% t %>% unlist
