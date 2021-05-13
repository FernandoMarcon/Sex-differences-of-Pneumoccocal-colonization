rm(list = ls())
pkgs <- c('tidyverse','rstatix','cluster','factoextra','ggrepel')
suppressPackageStartupMessages(sapply(pkgs, require, character.only = T))

#--- Settings
basedir <- 'intermediate/volunteer_wise_analysis'
gtm.db <- 'WikiPathways_2019_Human'

#--- Load Data
pheno <- read.delim(file.path(basedir, 'logFC_pheno.csv'))
head(pheno)

nes <- read.delim(file.path(basedir,'ssGSEA', gtm.db, 'NES_logFC_geneClean.csv'), row.names = 1)
colnames(nes) <- gsub('\\.','\\/',gsub('X','',colnames(nes)))
head(nes)

#--- K-means
nes.t <- nes %>% t %>% as.data.frame
fviz_nbclust(nes.t, kmeans, method='silhouette')
nes_kmeans <- kmeans(nes.t, centers = 2, nstart = 10)
fviz_cluster(nes_kmeans, data = nes.t)


pheno$cluster <- nes_kmeans$cluster[pheno$volunteer_id]
cluster.summary = pheno %>% separate(class, c('dataset','carriage','sex'), sep = '_') %>%
  unite('class', carriage, sex, sep = '_') %>% group_by(cluster,class) %>% summarize(num_vol = n())

ggplot(cluster.summary, aes(class, num_vol, fill = class)) + geom_bar(stat = 'identity') + facet_grid(.~as.factor(cluster)) +
  theme_linedraw() + theme(legend.position = 'top')

nes.full <- nes %>% rownames_to_column('pathway') %>% gather('volunteer_id','nes',-pathway)
nes.full = nes.full %>% merge(pheno,., by = 'volunteer_id',all.y = T)
nes.full = nes.full %>% separate(class, c('dataset','carriage','sex'), sep = '_')
nes.full = nes.full %>% separate(pathway, c('pathway','path_code'), sep = ' WP') %>% mutate(path_code = paste0('WP',path_code))
head(nes.full)

kurstal.test = nes.full %>% group_by(cluster, pathway) %>% kruskal_test(nes ~ sex) %>% mutate(.y. = gtm.db) %>% rename(DB = '.y.')
selected.pathways <- kurstal.test %>% filter(p < 0.01) %>% .$pathway %>% unique
plt.barplot <- kurstal.test %>% filter(pathway %in% selected.pathways) %>% mutate(p = -log10(p)) %>%
  ggplot(aes(p, reorder(pathway, p), fill = cluster)) + geom_bar(stat = 'identity', show.legend = F) + facet_grid(.~cluster) +
    labs(x = '-log10(p-value)',y = '', title = 'Female vs Male', subtitle = gtm.db) + theme_linedraw() +
    geom_vline(xintercept = -log10(.01), col= 'red', linetype = 'dashed') +
    geom_vline(xintercept = -log10(.05), col= 'red', linetype = 'dotted') +
    theme(panel.grid = element_blank())
plt.barplot

path.clusters <- kurstal.test %>% select(pathway, cluster, p) %>% filter(p < 0.1) %>% mutate(p = -log10(p)) %>%
  spread(cluster,p, fill = 1) %>% column_to_rownames('pathway') %>%
  setNames(paste0('cluster',colnames(.))) %>% rownames_to_column('pathway')

pdf(file.path(basedir, 'ssGSEA',gtm.db,'selectedPathways_byCluster_barplot.pdf'))
plt.barplot
# ggplot(path.clusters, aes(cluster1, cluster2)) + geom_point() + geom_text_repel(aes(label = pathway))
lapply(selected.pathways, function(pathway.name){
  nes.full %>% filter(pathway == pathway.name) %>% mutate(cluster = as.factor(cluster)) %>%
    ggplot(aes(sex,nes, fill = sex)) + geom_boxplot(show.legend =F, alpha = .5)   + geom_jitter(show.legend =F, aes(col = cluster)) +
      facet_grid(.~carriage) + theme_linedraw() + labs(title = pathway.name) + theme(panel.grid = element_blank())
  })
dev.off()

nes.full %>% filter(pathway == selected.pathways[1]) %>% mutate(cluster = paste0('cluster',cluster)) %>%
  spread(cluster,nes)

nes.full %>% filter(pathway %in% selected.pathways) %>% mutate(cluster = as.factor(cluster), group = paste0(carriage, '_',sex)) %>%
  ggplot(aes(group,nes, fill = group, col = group)) + geom_boxplot(show.legend =F)   + https://masternodes.online/geom_jitter(show.legend =F) +
    facet_grid(.~pathway) + theme_linedraw()



#--- Cross with luminex values
luminex <- read.delim(file.path(basedir, 'logFC_luminex.csv'))
head(luminex)

cytokine <- luminex %>% select(volunteer_id, timepoint, IL_7) %>%
  merge(pheno,.,by='volunteer_id',all.y=T)
head(cytokine)
ggplot(cytokine, aes(sex, IL_7)) + geom_boxplot() + geom_point() +
  facet_grid(.~carriage)

#### ===== BiClustering +++++ ####
# https://cran.r-project.org/web/packages/biclustermd/vignettes/Airports.html

# library(biclustermd)
# bc <- biclustermd(data = nes, col_clusters = 2, row_clusters = 4)
# autoplot(bc) +
#   scale_fill_gradient2(low = 'blue', high = 'red', mid = 'white') +
#   labs(x = 'Volunteers', y = 'Pathways', fill = 'NES')
#
#   # autoplot(bc$SSE)
#   # autoplot(bc$Similarities)
#
# logFC <- read.delim(file.path(basedir, 'logFC_geneClean.csv'), row.names = 1)
# colnames(logFC) <- gsub('\\.','\\/', gsub('X','',colnames(logFC)))
# head(logFC)
# bc <- biclustermd(data = logFC, col_clusters = 5, row_clusters = 4)
# autoplot(bc) +
#   scale_fill_gradient2(low = 'blue', high = 'red', mid = 'white') +
#   labs(x = 'Volunteers', y = 'Pathways', fill = 'logFC')
