rm(list = ls())
pkgs <- c('tidyverse')
suppressPackageStartupMessages(sapply(pkgs, require, character.only = T))
gtm.db = 'WikiPathways_2019_Human'
basedir <- file.path('intermediate/volunteer_wise_analysis/ssGSEA',gtm.db)

pathways <- c('IL−7 Signaling Pathway WP205',
              'Hematopoietic Stem Cell Gene Regulation by GABP alpha/beta Complex WP3657',
              'T−Cell antigen Receptor',
              'MicroRNAs in cardiomyocyte hypertrophy WP1544',
              'Mitochondrial Gene Expression WP391',
              'Interferon type I signaling pathways WP585',
              'Viral Acute Myocarditis WP4298',
              'T−Cell Receptor and Co−stimulatory Signaling WP2583',
              'Physiological and Pathological Hypertrophy of the Heart WP1528')
pathway = pathways[1]
# LOAD DATA
pheno <- read.delim('~/Documents/work/Sex-differences-of-Pneumoccocal-colonization/intermediate/volunteer_wise_analysis/logFC_pheno.csv')
head(pheno)
nes <- read.delim(file.path(basedir, 'NES_logFC_geneClean.csv'), row.names = 1)
colnames(nes) <- gsub('X','',gsub('\\.','\\/',colnames(nes)))
identical(colnames(nes), pheno$volunteer_id)
nes[pathway,] %>% t %>% as.data.frame
which(rownames(nes) == pathway)
pathway
rownames(nes)[grep('$WP205^', rownames(nes))]
