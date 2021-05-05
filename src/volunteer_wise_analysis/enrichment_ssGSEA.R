rm(list = ls())
pkgs <- c('tidyverse')
suppressPackageStartupMessages(sapply(pkgs, require, character.only = T))
# group.names <- c('POS_M','POS_F','NEG_M','NEG_F')
# dataset.names <- c('Adults1','Adults2','Adults3','Elderly1')
basedir = 'intermediate/volunteer_wise_analysis'
gtm.dbs <- list.files('data/GTM')

#### ============================================= ssGSEA ============================================= ####
gtm.db <- gtm.dbs[1]
outdir <- file.path(basedir, 'ssGSEA', gtm.db)
if(!dir.exists(outdir)) dir.create(outdir, recursive = F)

logFC.df = read.delim(file.path(basedir, 'logFC_data.csv'))
colnames(logFC.df) = gsub('X','',gsub('\\.','\\/',colnames(logFC.df)))

pheno <- read.delim(file.path(basedir, 'logFC_pheno.csv'))
rownames(pheno) <- pheno$volunteer_id

basedir <- '/home/marcon/Documents/work/Sex-differences-of-Pneumoccocal-colonization/'
setwd(file.path(basedir,'intermediate/volunteer_wise_analysis/ssGSEA/'))

# run Single_Sample_GSEA_ssGSEA_fgsea.R
source(file.path(basedir,'src/Single_Sample_GSEA_ssGSEA_fgsea.R'))
gmtfile <- file.path(basedir,'data/KEGG_2019_Human') # Reactome_2016
fileranks <- "logFC_ssGSEAinput.csv"
Ptype <- "padj"
pval_cutoff <- 0.1
ssGSEA(gmtfile=gmtfile,fileranks=fileranks,Ptype=Ptype,pval_cutoff=pval_cutoff)
