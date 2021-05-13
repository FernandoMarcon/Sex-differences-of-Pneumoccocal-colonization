rm(list = ls())
pkgs <- c('tidyverse', 'data.table')
suppressPackageStartupMessages(sapply(pkgs, require, character.only = T))


#--- Find all Luminex data
basedir <- '~/Documents/work/PhD/data/raw_data'
data <- read.delim(file.path(basedir, 'MIXED','Luminex.tsv'), row.names = 'sample_id') %>% select(-volunteer_id, -Day, -timepoint)
cytokine.names <- toupper(colnames(data))

all.filenames <- list.files(basedir, recursive = T, full.names = T)
all.filenames <- grep('tsv|csv|txt',all.filenames, value = T)
all.filenames <- all.filenames[-grep('count|rnaseq',all.filenames)]
all.filenames <- all.filenames[-grep('AGES_LAIV1_LAIV2_DATASHEET_DESSI_June_06_2019.csv',all.filenames)]

luminex.files <- lapply(all.filenames, function(filename){
  df <- fread(filename, nrows = 1)
  idx = sapply(cytokine.names, function(x) grep(paste0('^',x,'$'), colnames(df), ignore.case = T, value= T)) %>% unlist
  if(length(idx) == length(cytokine.names)) return(filename)
}) %>% unlist

#--- Combine all Luminex data
luminex.all <- lapply(luminex.files, function(luminex.fname) {
  data <- read.delim(luminex.fname, row.names = 'sample_id')
  idx = sapply(cytokine.names, function(x) grep(paste0('^',x,'$'), colnames(data), ignore.case = T))
  data <- data[,idx] %>% rownames_to_column('sample_id') %>% gather('cytokine', 'value', -sample_id, na.rm = T) %>%
    mutate(cytokine = toupper(cytokine))
  colnames(data)[3] <- basename(luminex.fname)
  data
  })
# lapply(luminex.all, head)

cytokine.all <- Reduce(function(x,y) merge(x, y, by = c('sample_id', 'cytokine'), all = T), luminex.all) %>%
  gather('table', 'value', -sample_id, -cytokine, na.rm = T) %>%
  group_by(sample_id, cytokine) %>% summarize(min_val = min(value)) %>% spread(cytokine,-sample_id)
write.table(cytokine.all, file.path('data','tidy_data','luminex.csv'), sep ='\t', quote = F, row.names = F)
