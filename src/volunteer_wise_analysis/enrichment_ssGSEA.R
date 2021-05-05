rm(list = ls())
pkgs <- c('tidyverse')
suppressPackageStartupMessages(sapply(pkgs, require, character.only = T))
basedir = 'intermediate/volunteer_wise_analysis'
gtm.dbs <- list.files('data/GTM')

#### ============================================= ssGSEA ============================================= ####
gtm.db <- gtm.dbs[1]
outdir <- file.path(basedir, 'ssGSEA', gtm.db)
if(!dir.exists(outdir)) dir.create(outdir, recursive = T)

file.copy(file.path(basedir, 'logFC_geneClean.csv'),outdir)
setwd(outdir)

# run Single_Sample_GSEA_ssGSEA_fgsea.R
prjt_path <- '/home/marcon/Documents/work/Sex-differences-of-Pneumoccocal-colonization'
source(file.path(prjt_path,'src/Single_Sample_GSEA_ssGSEA_fgsea.R'))
gmtfile <- file.path(prjt_path,'data/GTM/',gtm.db) # Reactome_2016
fileranks <- "logFC_geneClean.csv"
Ptype <- "padj"
pval_cutoff <- 0.1
ssGSEA(gmtfile=gmtfile,fileranks=fileranks,Ptype=Ptype,pval_cutoff=pval_cutoff)
file.remove('logFC_geneClean.csv')
