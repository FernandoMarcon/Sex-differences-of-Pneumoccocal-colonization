

#### =============== ssGSEA =============== ####
dir.create('intermediate/volunteer_wise_analysis/ssGSEA', showWarnings = F)

logFC.df = read.delim('intermediate/volunteer_wise_analysis/logFC.csv')
colnames(logFC.df) = gsub('X','',gsub('\\.','\\/',colnames(logFC.df)))

pheno <- read.delim('intermediate/volunteer_wise_analysis/logFC_pheno.csv')
rownames(pheno) <- pheno$volunteer_id

# ensembl to gene symbol
gene.dic <- read.delim('data/tidy_data/gene_annotation.csv')
temp = merge(gene.dic, logFC.df, by = 'gene_id')

# remove duplicated genes
temp = temp %>% filter(symbol != '') %>%
  mutate(gene_mean = rowMeans(.[,-c(1,2)])) %>%
  group_by(symbol) %>% filter(gene_mean == max(gene_mean)) %>%
  select(-gene_mean) %>% select(-gene_id) %>%
  rename(symbol = 'genes')
head(temp)
write.table(temp, 'intermediate/volunteer_wise_analysis/ssGSEA/logFC_ssGSEAinput.csv', sep = '\t',row.names = F, quote = F)

basedir <- '/home/marcon/Documents/work/Sex-differences-of-Pneumoccocal-colonization/'
setwd(file.path(basedir,'intermediate/volunteer_wise_analysis/ssGSEA/'))

# run Single_Sample_GSEA_ssGSEA_fgsea.R
source(file.path(basedir,'src/Single_Sample_GSEA_ssGSEA_fgsea.R'))
gmtfile <- file.path(basedir,'data/KEGG_2019_Human') # Reactome_2016
fileranks <- "logFC_ssGSEAinput.csv"
Ptype <- "padj"
pval_cutoff <- 0.1
ssGSEA(gmtfile=gmtfile,fileranks=fileranks,Ptype=Ptype,pval_cutoff=pval_cutoff)
