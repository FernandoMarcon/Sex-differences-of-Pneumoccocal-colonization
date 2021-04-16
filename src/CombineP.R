
df = pval.df

x = 4
values <- as.numeric(df[x,])
values[!is.na(values)]


sumlog(values)

library(metap)
df <- read.table(file = "C:/Users/hnaka/Desktop/Papers/Coronavirus_ACE/Mechanism/Pvalues.txt", 
                 header=TRUE, sep="\t", row.names = 1)

Pcombined <- unlist(lapply(rownames(df), function(x) { sumlog(as.numeric(df[x,])[!is.na(as.numeric(df[x,]))])$p } ))
df$FisherMethodP <- Pcombined

df$FDR <- p.adjust(Pcombined, method = "fdr", n = length(Pcombined))


write.table(df,file="FisherMethodPCombined.txt",sep="\t",col.names=NA)
                    
