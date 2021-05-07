rm(list = ls())
pkgs <- c('tidyverse')
suppressPackageStartupMessages(sapply(pkgs, require, character.only = T))

basedir <- 'intermediate/volunteer_wise_analysis'
gtm.db <- 'WikiPathways_2019_Human'

pheno <- read.delim(file.path(basedir, 'logFC_pheno.csv'))

nes <- read.delim(file.path(basedir,'ssGSEA', gtm.db, 'NES_logFC_geneClean.csv'), row.names = 1)
colnames(nes) <- gsub('\\.','\\/',gsub('X','',colnames(nes)))

#--- K-means
library(cluster)
library(factoextra)
nes.t <- nes %>% t %>% as.data.frame
nes_kmeans <- kmeans(nes.t, centers = 2, nstart = 25)
fviz_cluster(nes_kmeans, data = nes.t)

pheno$cluster <- nes_kmeans$cluster[pheno$volunteer_id]
cluster.summary = pheno %>% separate(class, c('dataset','carriage','sex'), sep = '_') %>%
  unite('class', carriage, sex, sep = '_') %>% group_by(cluster,class) %>% summarize(num_vol = n())
ggplot(cluster.summary, aes(class, num_vol, fill = class)) + geom_bar(stat = 'identity') + facet_grid(.~as.factor(cluster))

nes.full <- nes %>% rownames_to_column('pathway') %>% gather('volunteer_id','nes',-pathway)
nes.full = nes.full %>% merge(pheno,., by = 'volunteer_id',all.y = T)
nes.full = nes.full %>% separate(class, c('dataset','carriage','sex'), sep = '_')
head(nes.full)

library(rstatix)
kurstal.test = nes.full %>% group_by(cluster, pathway) %>% kruskal_test(nes ~ sex) %>% mutate(.y. = gtm.db) %>% rename(DB = '.y.')
selected.pathways <- kurstal.test %>% filter(p < 0.01) %>% .$pathway %>% unique
kurstal.test %>% filter(pathway %in% selected.pathways) %>% mutate(p = -log10(p)) %>%
  ggplot(aes(p, reorder(pathway, p), fill = cluster)) + geom_bar(stat = 'identity', show.legend = F) + facet_grid(.~cluster)
head(nes.full)

nes.full %>% filter(pathway == selected.pathways[1]) %>% mutate(cluster = as.factor(cluster)) %>%
  ggplot(aes(sex,nes)) + geom_boxplot(show.legend =F)   + geom_jitter(show.legend =F, aes(col = cluster)) +
    facet_grid(dataset~carriage) + theme_linedraw()

library(biclustermd)
bc <- biclustermd(data = nes, col_clusters = 2, row_clusters = 4)
autoplot(bc) +
  scale_fill_gradient2(low = 'blue', high = 'red', mid = 'white') +
  labs(x = 'Volunteers', y = 'Pathways', fill = 'NES')

  # autoplot(bc$SSE)
  # autoplot(bc$Similarities)

logFC <- read.delim(file.path(basedir, 'logFC_geneClean.csv'), row.names = 1)
colnames(logFC) <- gsub('\\.','\\/', gsub('X','',colnames(logFC)))
head(logFC)
bc <- biclustermd(data = logFC, col_clusters = 5, row_clusters = 4)
autoplot(bc) +
  scale_fill_gradient2(low = 'blue', high = 'red', mid = 'white') +
  labs(x = 'Volunteers', y = 'Pathways', fill = 'logFC')
