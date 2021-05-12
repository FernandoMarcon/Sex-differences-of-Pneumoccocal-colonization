rm(list = ls())
pkgs <- c('tidyverse', 'data.table')
suppressPackageStartupMessages(sapply(pkgs, require, character.only = T))


#--- Find all Luminex data
basedir <- '~/Documents/work/PhD/data/raw_data'
data <- read.delim(file.path(basedir, 'MIXED','Luminex.tsv'), row.names = 'sample_id') %>% select(-volunteer_id, -Day, -timepoint)
cytokine.names <- toupper(colnames(data))

all.filenames <- list.files(basedir, recursive = T, full.names = T)
all.filenames <- grep('tsv|csv|txt',all.filenames, value = T)

luminex.files <- lapply(all.filenames, function(filename){
  df <- fread(filename, nrows = 1)
  idx = sapply(cytokine.names, function(x) grep(paste0('^',x,'$'), colnames(df), ignore.case = T))
  if(length(idx) == length(cytokine.names)) return(filename)
}) %>% unlist
length(luminex.files)
sapply(basename(luminex.files), print)

#--- Combine all Luminex data
luminex.fname <- luminex.files[1]
luminex.all <- lapply(luminex.files, function(luminex.fname) {
  message(luminex.fname)
  data <- read.delim(luminex.fname, row.names = 'sample_id')
  idx = sapply(cytokine.names, function(x) grep(paste0('^',x,'$'), colnames(data), ignore.case = T))
  data[,idx]
  })
head()
colnames(data) <- toupper(colnames(data))





logFC <- read.delim('intermediate/volunteer_wise_analysis/logFC_pheno.csv')
head(logFC)

intersect(logFC$volunteer_id, unique(data$volunteer_id))
setdiff(logFC$volunteer_id, unique(data$volunteer_id))
setdiff(unique(data$volunteer_id), logFC$volunteer_id)
unique(data$volunteer_id)



laiv2 <- read.delim(file.path(basedir, 'LAIV2','luminex', 'LAIV2_luminex_data.tsv'), row.names = 'sample_id')
colnames(laiv2)
pilot_laiv1 <- read.delim(file.path(basedir, 'LAIV1', 'LAIV1_Luminex_values.tsv'), row.names = 'sample_id') %>%
  select(-volunteer_id, -study, -day)
head(pilot_laiv1)

cytokines <- colnames(laiv2)
laiv2 <- laiv2[,cytokines]
pilot_laiv1 <- pilot_laiv1[,cytokines]

cytokines.df <- rbind(laiv2, pilot_laiv1)
head(cytokines.df)

volunteer_ids <- gsub('_.*', '', rownames(cytokines.df))
setdiff(logFC$volunteer_id, volunteer_ids)
