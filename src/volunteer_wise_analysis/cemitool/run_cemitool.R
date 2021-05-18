rm(list = ls())
pkgs <- c('tidyverse', 'CEMiTool')
suppressPackageStartupMessages(sapply(pkgs, require, character.only = T))

basedir <- 'intermediate/volunteer_wise_analysis'
outdir <- file.path(basedir, 'cemitool')

data <- read.delim(file.path(basedir, 'logFC_geneClean.csv'), row.names = 'genes')
colnames(data) <- gsub('\\.','\\/',gsub('X','',colnames(data)))
data[1:4,1:4]

pheno <- read.delim(file.path(basedir, 'logFC_pheno.csv'))
head(pheno)

# GMT
gmt <- read_gmt('data/GTM/WikiPathways_2019_Human')

# Interactions
int_fname <- system.file("extdata", "interactions.tsv", package = "CEMiTool")
int_df <- read.delim(int_fname)
head(int_df)

#---
comp = 'dataset_carriage_sex'

cem <- cemitool(expr = data, annot = pheno, sample_name_column="volunteer_id", class_column="class",
                gmt = gmt, interactions = int_df, filter = T, plot = T, verbose = T)

generate_report(cem, directory=file.path(outdir, comp, "./Report"))
write_files(cem, directory=file.path(outdir, comp, "./Tables"))
save_plots(cem, "all", directory=file.path(outdir, comp, "./Plots"))

plot_gsea(cem)
gsea.plt <- show_plot(cem, "gsea")[[1]] + theme(axis.text.x = element_text(angle = 90))
pdf(file.path(outdir, comp, 'Plots','gsea_fixed.pdf'))
gsea.plt
dev.off()

#---
comp = 'carriage_sex'

temp = pheno %>% separate(class, c('dataset', 'carriage', 'sex'), sep = '_') %>%
    unite('class', carriage, sex)

cem <- cemitool(expr = data, annot = temp, sample_name_column="volunteer_id", class_column="class",
                gmt = gmt, interactions = int_df, filter = T, plot = T, verbose = T)

generate_report(cem, directory=file.path(outdir, comp, "./Report"))
write_files(cem, directory=file.path(outdir, comp, "./Tables"))
save_plots(cem, "all", directory=file.path(outdir, comp, "./Plots"))

gsea.plt <- show_plot(cem, "gsea")[[1]] + theme(axis.text.x = element_text(angle = 90))
pdf(file.path(outdir, comp, 'Plots','gsea_fixed.pdf'))
gsea.plt
dev.off()

#---
comp = 'sex'

temp = pheno %>% separate(class, c('dataset', 'carriage', 'sex'), sep = '_')

cem <- cemitool(expr = data, annot = temp, sample_name_column="volunteer_id", class_column="sex",
                gmt = gmt, interactions = int_df, filter = T, plot = T, verbose = T)

generate_report(cem, directory=file.path(outdir, comp, "./Report"))
write_files(cem, directory=file.path(outdir, comp, "./Tables"))
save_plots(cem, "all", directory=file.path(outdir, comp, "./Plots"))

gsea.plt <- show_plot(cem, "gsea")[[1]] + theme(axis.text.x = element_text(angle = 90))
pdf(file.path(outdir, comp, 'Plots','gsea_fixed.pdf'))
gsea.plt
dev.off()
