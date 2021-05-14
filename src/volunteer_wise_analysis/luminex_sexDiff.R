rm(list = ls())
pkgs <- c('tidyverse','rstatix','cluster','factoextra','ggrepel')
suppressPackageStartupMessages(sapply(pkgs, require, character.only = T))

#--- SETTINGS
basedir <- 'intermediate/volunteer_wise_analysis'

#--- DATA
luminex <- read.delim(file.path(basedir, 'logFC_luminex.csv'))
pheno <-  read.delim(file.path(basedir, 'logFC_pheno.csv'))

#--- TIDY
pheno <- merge(luminex[,1:2],pheno, by = 'volunteer_id')
luminex <- luminex %>% unite('sample_id', volunteer_id, timepoint, sep = '_') %>% column_to_rownames('sample_id')
pheno <- pheno %>% unite('sample_id', volunteer_id, timepoint, sep = '_', remove = F) %>% column_to_rownames('sample_id') %>%
  separate(class, c('dataset','carriage','sex'), sep = '_')
identical(rownames(pheno), rownames(luminex))

head(pheno)

table(pheno$timepoint)
kurstal.test = cbind(pheno[,c('timepoint','sex')], luminex) %>%  gather('cytokine','value', -timepoint, -sex) %>%
  group_by(timepoint, cytokine) %>% kruskal_test(value ~ sex)
head(kurstal.test)

kurstal.test %>% filter(p < 0.1) %>%
  ggplot(aes(-log10(p), cytokine, fill = cytokine)) + geom_bar(stat = 'identity', show.legend = F) +
  facet_grid(timepoint~., scales = 'free') + theme_linedraw(15) +
  geom_vline(xintercept = -c(log10(.1),log10(.05),log10(.01)), linetype = 'dashed', col = 'red') +
  theme(panel.grid = element_blank())


head(pheno)
kurstal.test = cbind(pheno[,c('dataset','timepoint','sex')], luminex) %>%  gather('cytokine','value', -dataset, -timepoint, -sex) %>%
  group_by( dataset, timepoint, cytokine) %>% kruskal_test(value ~ sex)
head(kurstal.test)

kurstal.test %>% filter(p < .1) %>%
  ggplot(aes(-log10(p), cytokine, fill = cytokine)) + geom_bar(stat = 'identity', show.legend = F) +
  facet_grid(timepoint~dataset, scales = 'free', space = 'free') + theme_linedraw(15) +
  geom_vline(xintercept = -c(log10(.1),log10(.05),log10(.01)), linetype = 'dashed', col = 'red') +
  theme(panel.grid = element_blank())

luminex %>% gather('cytokine', 'value', na.rm = T) %>% ggplot(aes(value, group = cytokine, col = cytokine)) + geom_density(show.legend = F)

prcomp(luminex)$x %>% data.frame %>% ggplot(aes(PC1,PC2)) + geom_point()
