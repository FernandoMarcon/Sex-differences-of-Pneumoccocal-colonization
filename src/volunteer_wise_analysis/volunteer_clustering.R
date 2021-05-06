rm(list = ls())
pkgs <- c('tidyverse')
suppressPackageStartupMessages(sapply(pkgs, require, character.only = T))

basedir <- 'intermediate/volunteer_wise_analysis'
gtm.db <- 'WikiPathways_2019_Human'

nes <- read.delim(file.path(basedir,'ssGSEA', gtm.db, 'NES_logFC_geneClean.csv'), row.names = 1)
colnames(nes) <- gsub('\\.','\\/',gsub('X','',colnames(nes)))
head(nes)

library(cluster)
library(factoextra)
nes.t <- nes %>% t %>% as.data.frame
nes_kmeans <- kmeans(nes.t, centers = 3, nstart = 25)
fviz_cluster(nes_kmeans, data = nes.t)


library(biclustermd)
bc <- biclustermd(data = nes, col_clusters = 5, row_clusters = 4)
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
