rm(list = ls())
pkgs <- c('tidyverse','pheatmap','RColorBrewer','rstatix')
suppressPackageStartupMessages(sapply(pkgs, require, character.only = T))

#### =============== FUNCTIONS ================ ####
loadEnrichementData <- function(gtm.db, basedir, p.thrs = 0.01, vol.perc = .7) {
  nes <- read.delim(file.path(basedir, 'NES_logFC_geneClean.csv'), row.names = 1)
  padj <- read.delim(file.path(basedir, 'padj_logFC_geneClean.csv'), row.names = 1)
  colnames(padj) <- colnames(nes) <- gsub('X','',gsub('\\.','\\/',colnames(nes)))

  #--- FILTERING
  # Filter by padj < threshold for at least X% of volunteers in a given group
  selected.pathways <- lapply(group.names, function(group.name){ #group.name = group.names[1]
    selected.vol = pheno %>% filter(group == group.name) %>% rownames
    temp <- padj[,intersect(colnames(padj), selected.vol)]
    selected.pathways <- names(rowSums(temp < p.thrs) > ncol(temp)*vol.perc)
    selected.pathways
    }) %>% unlist %>% unique

  return(nes[selected.pathways, ])
    # return(nes)
}

# dataset.names <- c('Adults1','Adults2','Adults3','Elderly1')
group.names <- c('POS_M','POS_F','NEG_M','NEG_F')
gtm.dbs <- list.files('data/GTM')

#--- LOAD DATA
pheno <- read.delim(file.path('intermediate/volunteer_wise_analysis/logFC_pheno.csv'), row.names = 'volunteer_id') %>%
  separate(class, c('dataset', 'carriage','sex'), sep = '_') %>%
  unite('group', carriage, sex, sep = '_', remove = F)

#--- Kruskal-Wallis rank sum test
krustal.all <- lapply(gtm.dbs, function(gtm.db) {#  gtm.db = gtm.dbs[1]
  basedir <- file.path('intermediate/volunteer_wise_analysis/ssGSEA',gtm.db)

  nes <- loadEnrichementData(gtm.db, basedir)

  nes.full = nes %>% rownames_to_column('pathway') %>% gather('volunteer_id','nes',-pathway)
  nes.full = pheno %>% rownames_to_column('volunteer_id') %>% merge(., nes.full, by = 'volunteer_id', all.y = T)
  nes.full = nes.full %>% group_by(dataset, pathway) %>% kruskal_test(nes ~ sex) %>% mutate(.y. = gtm.db) %>% rename(DB = '.y.')
  # write.table(nes.full, file.path(basedir, 'sexDiff_kurstallTest.csv'), sep ='\t', row.names = F, quote =F)
  nes.full
})

pval.thrs = 0.1
temp = Reduce(rbind, krustal.all) %>% filter(p < pval.thrs) %>% mutate(p = -log10(p)) %>%
        mutate(pathway = gsub(' \\(.*','',pathway), pathway = gsub(' Homo.*','',pathway),
                        DB = gsub('GO_Biological_Process_2018','GO-BP',DB),
                        DB = gsub('GO_Cellular_Component_2018','GO-CC',DB),
                        DB = gsub('GO_Molecular_Function_2018','GO-MF',DB),
                        DB = gsub('KEGG_2019_Human','KEGG',DB),
                        DB = gsub('Reactome_2016','Reactome',DB),
                        DB = gsub('WikiPathways_2019_Human','WikiPathways',DB),
                      dataset = gsub('Adults','A',dataset))
plt.krustal <- ggplot(temp, aes(p,reorder(pathway, p), fill = DB)) + geom_bar(stat = 'identity', show.legend = F) +
  facet_grid(DB~dataset,scales = 'free', space = 'free_y') +
  labs(y = '', x = '-log10(pvalue)', title = 'Female vs Male', subtitle = 'Krustal-Wallis Test') +
  theme_linedraw() + theme(panel.grid = element_blank(),axis.text=element_text(size=2),axis.title=element_text(size=2))

pdf(paste0('intermediate/volunteer_wise_analysis/ssGSEA/sexDiff_kustalTest_selectedPathways_p',pval.thrs,'.pdf'),height = 40, width = 4)
plt.krustal
dev.off()

# SELECT DB
lapply(gtm.dbs, function(gtm.db) {
  gtm.db = gtm.dbs[1]
  basedir <- file.path('intermediate/volunteer_wise_analysis/ssGSEA',gtm.db)


  plt.nes <- pheatmap(nes, show_rownames = F, show_colnames = F, annotation_col = pheno,
    main = paste0('Pathways with padj < ',p.thrs,' for > ',vol.perc*100,'% of volunteers in a group'))

  pdf(file.path(basedir, paste0('NES_selectedPathways_padj_',p.thrs,'_volPerc_',vol.perc*100,'_heatmap.pdf')))
  plt.nes
  dev.off()


  #--- T-test: Female vs Male
  female.vols <- pheno %>% filter(sex == 'F') %>% rownames
  male.vols <- pheno %>% filter(sex == 'M') %>% rownames

  tTest_res <- apply(nes, 1, function(x) t.test(abs(x[female.vols]), abs(x[male.vols]))$p.value) %>%
    as.data.frame %>% setNames('pval') %>% rownames_to_column('pathway')

  write.table(tTest_res, file.path(basedir, paste0(gtm.db, '_selectedPathways.csv')), sep = '\t', quote =  F, row.names = F)

  plt.pathways <- tTest_res %>% filter(pval < 0.01) %>% mutate(pval = -log10(pval), pathway = gsub(' \\(.*','',pathway)) %>%
    ggplot(aes(pval, reorder(pathway,pval), fill = pval)) + geom_bar(show.legend = F, stat = 'identity') +
      labs(x = '-log10(Pvalue)', y = '', title = gtm.db, subtitle = 'pathways with pval < 0.01')

  pdf(file.path(basedir, paste0(gtm.db, '_selectedPathways_barplot.pdf')),width = 10)
  print(plt.pathways)
  dev.off()
  })
