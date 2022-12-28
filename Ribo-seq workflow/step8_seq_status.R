#!/usr/bin/env Rscript

## take the arguments from environment
args = commandArgs(trailingOnly=TRUE)


# test if there is at least one argument: if not, return an error
if (length(args)!=1) {
  stop("Please supply the sepecies for the script \n At least one argument must be supplied (input file).n", call.=FALSE)}

#assign the args to the variates 
species =  args[1]

# use to generate graph to show the periodicty of sequencing, as well as the distribution of reads. 
library(dplyr)
library(reshape2)
library(ggplot2)
#================================== periodicity ==========================

#setting for test


# sam is sample name, up_st and down_st are the number of aa around start codon,
# end the number of AA you want you check around stop codon.
plot.per = function (sam, up_st, down_st, up_end, down_end)     {
  
  
  df = read.table(paste0(sam, "_A.bg"), sep = "\t")
  df = df[df$V4 != 0,]
  
  
  df1 = df %>% group_by(V1) %>% mutate(avg = V4/mean(V4))
  
  # check -29 to 29 around start codon
  
  dfs = df1[df1$V2 > (100 - up_st*3) & df1$V2 < (100 + down_st*3),]
  
  dfs1 = as.data.frame(aggregate(avg~V2, dfs, sum))
  dfs1$V2 = dfs1$V2 - 100
  dfs1$col = "white"
  dfs1$col[dfs1$V2 > 0 & dfs1$V2%%3 == 0] = "grey88"
  dfs1$col[dfs1$V2 > 0 & dfs1$V2%%3 == 1] = "lightpink2"
  dfs1$col[dfs1$V2 > 0 & dfs1$V2%%3 == 2] = "seagreen4"
  setEPS()
  #postscript(paste("codon_stat/", sam, "start periodicty.eps"), width = 10, height = 5)
  
  png(filename = paste("codon_stat/", sam, "start periodicty.png"), width = 950, height = 500, bg = "white")
  
  par(mar=c(5,6,5,4),mgp=c(4,1,0))
  
  barplot(dfs1$avg, names.arg = as.character(dfs1$V2), main = "start codon", col = dfs1$col,
          ylab = "Reads Count",xlab = "Position",cex.axis=2,cex.names = 2.5,cex.lab = 2.5,cex.main=3)
  
  dev.off()
  # check stop codon 
  
  ref = read.table(paste0("/media/hp/disk1/song/Genomes/", species,"/",species, "_ref/longest_cDNA.chrom"), sep = "\t")
  names(ref) = c("gene", "length")
  df1 = merge(ref, df1, by.x = "gene", by.y = "V1")
  df1$pos = df1$V2 - df1$length + 100
  
  dfe = df1[df1$pos < down_end*3 & df1$pos >= -up_end*3, ]
  
  dfe1 = as.data.frame(aggregate(avg~pos, dfe, sum))
  
  dfe1$col = "white"
  dfe1$col[dfe1$pos < 0 & dfe1$pos%%3 == 0] = "grey88"
  dfe1$col[dfe1$pos < 0 & dfe1$pos%%3 == 1] = "lightpink2"
  dfe1$col[dfe1$pos < 0 & dfe1$pos%%3 == 2] = "seagreen4"
  
  setEPS()
  #postscript(paste("codon_stat/", sam, "stop periodicty.eps"), width = 10, height = 5)
  png(filename = paste("codon_stat/", sam, "stop periodicty.png"), width = 950, height = 500, bg = "white")
  
  par(mar=c(5,6,5,4),mgp=c(4,1,0))
  
  barplot(dfe1$avg, names.arg = as.character(dfe1$pos), main = "stop codon", 
          col = dfe1$col,ylab = "Reads Count",xlab = "Position",cex.axis=2,cex.names = 2.5,cex.lab = 2.5,cex.main=3)
  
  dev.off()
  
}

line.density= function (type, sam)    {
  df = read.table(paste0("coverage_stat/", sam, "_", type, "_distribution_meta.txt" ), sep = "\t")

  m = as.matrix(df[,2:ncol(df)])
  rownames(m) = df[,1]
  d1 = as.data.frame( m/rowSums(m) )
  
  colnames(d1)[1:ncol(d1)] = -10:(ncol(d1)-10)
  d1$gene = rownames(d1)
  
  d2 = melt(d1, value.name = "Value", id.vars = "gene")
  
  d3 = d2 %>% group_by(variable) %>% summarise(mean = mean(Value), std = sd(Value) )
  
  eb = aes(ymax = mean + std, ymin = mean - std)
  p = ggplot(data = d3, aes(x=variable, y = mean, group = 1)) +
    geom_line(size = 0.5, color = "blue") +
    theme_bw() + 
    theme(axis.line = element_line(colour = "black"),
          plot.title = element_text(color="red", size=9, face="bold.italic"), # title 
          axis.title.x = element_text(color="blue", size=7, face="bold"),
          axis.title.y = element_text(color="#993333", size=7, face="bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          #axis.ticks = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) + 
    scale_x_discrete(breaks = seq(0,120,10),
                     )
  p = p + labs(y = "reads density", x= "position of gene" ) + ggtitle(paste(sam, type))
  
  
    #+ 
    #geom_ribbon(eb, alpha = 0.5) 
  ggsave(filename = paste0("coverage_stat/", sam, "_", type, "_meta.jpg"), plot=p, width = 10, height = 5, dpi = 200, units = "cm")
}




sample = list.files(path = "./", pattern = ".bg")
sample = gsub("_A.bg", "", sample)
#sample = gsub(".fq.gz", "", sample)

for (s in sample)  {
  #plot.per("RNA", s, 30, 30, 20, 20)
  plot.per( s, 20, 20, 20, 20)
  
  line.density("RPF", s)
  #line.density("RNA", s)
}

