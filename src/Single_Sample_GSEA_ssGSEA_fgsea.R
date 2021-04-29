################################################################
#Packages required: tidyr, dplyr, data.table, clusterProfiler
#                   fgsea,ggplot2
library(tidyr)
library(dplyr)
options(stringsAsFactors = F)
pkgs <- c('data.table', 'clusterProfiler', 'fgsea', 'ggplot2')
functioning <- sapply(pkgs, function(x) require(x, character.only=T))
if(!all(functioning)){stop('Install Required Packages')}

################################################################
################################################################
####Do not touch here. Just load functions######################
################################################################
################################################################
################################################################
read.gmt <- function(fname){
  res <- list(genes=list(), 
              desc=list())
  gmt <- file(fname)
  gmt.lines <- readLines(gmt)
  close(gmt)
  gmt.list <- lapply(gmt.lines, 
                     function(x) unlist(strsplit(x, split="\t")))
  gmt.names <- sapply(gmt.list, '[', 1)
  gmt.desc <- lapply(gmt.list, '[', 2)
  gmt.genes <- lapply(gmt.list,
                      function(x){x[3:length(x)]})
  names(gmt.desc) <- names(gmt.genes) <- gmt.names
  return(gmt.genes)
}

ssGSEA <- function(gmtfile=gmtfile,fileranks=fileranks,Ptype=Ptype,pval_cutoff=pval_cutoff){
  ####Do not touch. Just run#######################
  read.gmt <- function(fname){
    res <- list(genes=list(), 
                desc=list())
    gmt <- file(fname)
    gmt.lines <- readLines(gmt)
    close(gmt)
    gmt.list <- lapply(gmt.lines, 
                       function(x) unlist(strsplit(x, split="\t")))
    gmt.names <- sapply(gmt.list, '[', 1)
    gmt.desc <- lapply(gmt.list, '[', 2)
    gmt.genes <- lapply(gmt.list,
                        function(x){x[3:length(x)]})
    names(gmt.desc) <- names(gmt.genes) <- gmt.names
    return(gmt.genes)
  }
  
  #Load gene sets
  gmt <- read.gmt(gmtfile)
  #Load ranks
  ranks_df  <- read.delim(fileranks)
  #remove rows with NA in any column
  ranks_df <- ranks_df[complete.cases(ranks_df), ]
  colnames(ranks_df)[1] <- "Genes"
  
  ranks_ch <- colnames(ranks_df)[-1]
  ranks_ch <- setNames(ranks_ch, ranks_ch)
  
  #run fastGSEA
  tmp_ranks <- lapply(ranks_ch, function(rankname){
    tmpdf <- ranks_df[,c('Genes', rankname)]
    tmpdf <- tmpdf[complete.cases(tmpdf),]
    tmpranks <- tmpdf[[rankname]]
    names(tmpranks) <- tmpdf$Genes
    fgseaRes <- fgsea(pathways = gmt, stats = tmpranks, 
                      minSize = 15, maxSize = 2000, nperm = 1000)
    fgseaRes <- as.data.frame(fgseaRes)
    fgseaRes <- fgseaRes[,c('pathway', 'pval', 'padj', 'NES','size','leadingEdge')]
    fgseaRes
  })
  
  
  # Remove ranks without enrichment
  tmp_ranks <- Filter(function(x) nrow(x) > 1, tmp_ranks)
  
  # Write output NES
  rm(df)
  for(name_rank in names(tmp_ranks)){
    rank_out <- tmp_ranks[[name_rank]] 
    xxx   <- rank_out[,"NES"]
    names(xxx) <- rank_out[,"pathway"]
    df <- cbind(df,xxx)
    colnames(df)[ncol(df)] <- name_rank
    #rank_nameout <- paste0('fgsea_', name_rank, '.tsv')
    #write.table(rank_out, file = rank_nameout, sep = '\t', quote = F, row.names = F)
  }
  df <- df[,-1]
  rank_nameout <- paste0('NES_', fileranks,sep="")
  write.table(df, file = rank_nameout, sep = '\t', quote = F, row.names = T,col.names=NA)
  
  # Write output AdjP
  rm(df)
  for(name_rank in names(tmp_ranks)){
    rank_out <- tmp_ranks[[name_rank]] 
    xxx   <- rank_out[,"padj"]
    names(xxx) <- rank_out[,"pathway"]
    df <- cbind(df,xxx)
    colnames(df)[ncol(df)] <- name_rank
    #rank_nameout <- paste0('fgsea_', name_rank, '.tsv')
    #write.table(rank_out, file = rank_nameout, sep = '\t', quote = F, row.names = F)
  }
  df <- df[,-1]
  rank_nameout <- paste0('padj_', fileranks,sep="")
  write.table(df, file = rank_nameout, sep = '\t', quote = F, row.names = T,col.names=NA)
  
  # Write output pval
  rm(df)
  for(name_rank in names(tmp_ranks)){
    rank_out <- tmp_ranks[[name_rank]] 
    xxx   <- rank_out[,"pval"]
    names(xxx) <- rank_out[,"pathway"]
    df <- cbind(df,xxx)
    colnames(df)[ncol(df)] <- name_rank
    #rank_nameout <- paste0('fgsea_', name_rank, '.tsv')
    #write.table(rank_out, file = rank_nameout, sep = '\t', quote = F, row.names = F)
  }
  df <- df[,-1]
  rank_nameout <- paste0('pval_', fileranks,sep="")
  write.table(df, file = rank_nameout, sep = '\t', quote = F, row.names = T,col.names=NA)
  
  # Write output Leading Edge genes
  rm(df)
  for(name_rank in names(tmp_ranks)){
    rank_out <- tmp_ranks[[name_rank]] 
    xxx   <- rank_out[,"leadingEdge"]
    names(xxx) <- rank_out[,"pathway"]
    df <- cbind(df,xxx)
    colnames(df)[ncol(df)] <- name_rank
    #rank_nameout <- paste0('fgsea_', name_rank, '.tsv')
    #write.table(rank_out, file = rank_nameout, sep = '\t', quote = F, row.names = F)
  }
  df <- df[,-1]
  rank_nameout <- paste0('LE_', fileranks,sep="")
  write.table(df, file = rank_nameout, sep = '\t', quote = F, row.names = T,col.names=NA)
  
  nes <- data.table::fread(paste0('NES_', fileranks,sep=""))
  pval <- data.table::fread(paste0(Ptype,'_', fileranks,sep=""))
  
  nes_melt <- nes %>%
    rename(pathway=V1) %>%
    gather(sample, nes, -pathway)
  
  pval_melt <- pval %>%
    rename(pathway=V1) %>%
    gather(sample, pval, -pathway)
  
  result <- full_join(nes_melt, pval_melt, by=c("pathway", "sample")) %>%
    filter(pval <= pval_cutoff) %>%
    select(pathway, sample, nes) %>%
    spread(sample, nes)
  
  rank_nameout <- paste0('NES_',Ptype,pval_cutoff,"_", fileranks,sep="")
  write.table(result, file = rank_nameout, sep = '\t', quote = F, row.names = T,col.names=NA)
}
##################################################
##################################################
##################################################



################################################################
#Input files:
#Rank file = contains the Xv0 of each patient
#Gene set = Pathways in gmt format

################################################################
##################Change things here only#######################
################################################################
setwd("~/Downloads/")
#Gene set file name:
#(make sure there is no duplicated gene set name)
gmtfile <- "GeneSet_test.gmt"
#Rank file name:
fileranks <- "Xv0_test.txt"
#Type of P-value significance for NES enrichment
#Change Ptype to "pval" in case you want to use nominal P
#instead of Adjusted P-value
Ptype <- "padj"
#Cutoff used to call a pathway significant
pval_cutoff <- 0.1

#Run ssGSEA
ssGSEA(gmtfile=gmtfile,fileranks=fileranks,Ptype=Ptype,pval_cutoff=pval_cutoff)

