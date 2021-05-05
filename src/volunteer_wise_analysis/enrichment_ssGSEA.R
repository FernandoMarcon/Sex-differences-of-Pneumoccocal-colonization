rm(list = ls())
pkgs <- c('tidyverse')
suppressPackageStartupMessages(sapply(pkgs, require, character.only = T))

#### ============================================= ssGSEA ============================================= ####
prjt_path <- '/home/marcon/Documents/work/Sex-differences-of-Pneumoccocal-colonization'
basedir = file.path(prjt_path, 'intermediate','volunteer_wise_analysis')
gtm.dbs <- list.files(file.path(prjt_path, 'data','GTM'), full.names = T)

fileranks <- "logFC_geneClean.csv"
Ptype <- "padj"
pval_cutoff <- 0.1

# run Single_Sample_GSEA_ssGSEA_fgsea.R
lapply(gtm.dbs, function(gtm.db) {
  # gtm.db <- gtm.dbs[1]
  outdir <- file.path(basedir, 'ssGSEA', basename(gtm.db))
  if(!dir.exists(outdir)) dir.create(outdir, recursive = T)

  file.copy(file.path(basedir, 'logFC_geneClean.csv'), outdir)
  setwd(outdir)
  source(file.path(prjt_path, 'src/Single_Sample_GSEA_ssGSEA_fgsea.R'))

  ssGSEA(gmtfile=gtm.db,fileranks=fileranks,Ptype=Ptype,pval_cutoff=pval_cutoff)
  file.remove('logFC_geneClean.csv')

  })
