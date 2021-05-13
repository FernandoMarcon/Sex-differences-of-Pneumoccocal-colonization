rm(list = ls())
pkgs <- c('tidyverse')
suppressPackageStartupMessages(sapply(pkgs, require, character.only = T))
basedir <- 'intermediate/volunteer_wise_analysis'

pheno <- read.delim(file.path(basedir,'logFC_pheno.csv'))
head(pheno)
luminex <- read.delim('data/tidy_data/luminex.csv') %>%
  separate(sample_id, c('volunteer_id', 'timepoint'), sep = '_') %>%
  filter(volunteer_id %in% unique(pheno$volunteer_id))
luminex[,-c(1:2)] <- log10(luminex[,-c(1:2)] + .5)
head(luminex)

# temp = luminex[,1:2] %>% mutate(study = sapply(strsplit(volunteer_id, '\\/'), `[[`, 2))
# temp %>% group_by(study, timepoint) %>% summarize(n())

boxplot(luminex[,-c(1:2)])


luminex %>% gather(cytokine, value, -volunteer_id, -timepoint)

luminex.logFC <- luminex %>% gather(cytokine, value, -volunteer_id, -timepoint) %>%
  group_by(volunteer_id, cytokine) %>% mutate(bl = value[which(timepoint == 'baseline')]) %>%
  filter(complete.cases(value)) %>% mutate(logFC = value - bl) %>%
  filter(timepoint != 'baseline') %>% select(-value, -bl) %>%
  spread(cytokine, -volunteer_id)

write.table(luminex.logFC, file.path(basedir, 'logFC_luminex.csv'), sep = '\t', quote = F, row.names = F)
