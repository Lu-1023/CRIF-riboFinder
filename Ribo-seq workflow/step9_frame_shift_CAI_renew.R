#set the plot_sort for ploting CAI_CBI correlation or decode rate corrlation frameshift barplot
plot_sort="CAI_CBI"

if (plot_sort=="CAI_CBI"){
  setwd("~/Data/ribo_seq/chlamydomonas/coverage_stat")
}else{
  setwd("~/Data/ribo_seq/chlamydomonas/coverage_stat/RPF")
}

library(ggplot2)
samples<-system("ls ~/Data/ribo_seq/chlamydomonas/coverage_stat/*RPF_count.txt",intern = T)
for ( s in samples )    {
  plot.graph(s,plot_sort)
}


#read data 
plot.graph <- function(sam,sorts)  {
  
  d <- read.table(sam, sep="\t", header = T, stringsAsFactors = F)
  
  # calculate the number codon that are covered by at least 1 read 
  d$inframe_codon_cover_number = ceiling(d$inframe_cov/100 * (d$ORF_length - 200)/3 )
  # note  hist(d$inframe_codon_cover_number)
  d$RPF_total = apply(d[,c("inframe_count", "shift_plus_count", "shift_minus_count")], 1, sum)
  d$RPF_RPKM = d$RPF_total*1000000*1000/sum(d$RPF_total)/(d$ORF_length-200) 

  # note!! the y value
  if (sorts=="CAI_CBI"){
    d2 <- read.table("CBI_CAI_0.txt", sep = "\t", header=T, fill=T, stringsAsFactors  = F)
  }else{
    y<-gsub(".*_(.*)","RNASeq_\\1",sam)
    d2 <- read.table(paste0("decode_",y,"_100.txt"), sep = "\t", header=T, fill=T, stringsAsFactors  = F)
  }
  
  
  
  dfe <- merge(d2, d, by.x="transcription_id", by.y="Gene") 
  
  # dfe remove NA------------------------------------------------------------------------------------------------

  
  dfe$shift_plus_ratio = dfe$shift_plus_count/dfe$RPF_total
  dfe$shift_minus_ratio = dfe$shift_minus_count/dfe$RPF_total
  dfe$shift_ratio = (dfe$shift_minus_count + dfe$shift_plus_count)/dfe$RPF_total
  
  df <- dfe[dfe$inframe_cov > 5 & dfe$inframe_codon_cover_number > 10,]
  df<-transform(df,RNA_seq_cov=NULL)
  df<-transform(df,RNA_seq_count=NULL)
  df <- na.omit(df)
  
  df1.s <- dfe[order(dfe$shift_ratio),]
  df1.s<-df1.s[df1.s$RPF_RPKM > 30,]
  bin = 10
  
  for (i in 1:bin)  {
    range = as.integer(length(df1.s$shift_ratio)/bin)
    start = range*(i-1)+1
    end = range*i
    for (j in start:end)    {
      df1.s[j,"rank"] = i 
    }
    
  }
  df1.s$rank = as.numeric(df1.s$rank)
  
  
  # draw boxplot 
  jpeg(filename = paste(sam, " boxplot.jpeg", sep=""),
       width = 350, height = 450, units = "px", quality = 100)
  
  par(mfrow=c(1,1))
  if (sorts=="CAI_CBI"){
      boxplot(avg_CAI~rank, data=df1.s, col = rainbow(bin), axes = F, notch = T,
          main = paste(sam, "gene number =", length(df1.s$shift_ratio), sep = " "),
          xlab = "median of out/in ratio", ylab = "CAI")
  }else{
    boxplot(avg_deco~rank, data=df1.s, col = rainbow(bin), axes = F, notch = T,
          main = paste(sam, "gene number =", length(df1.s$shift_ratio), sep = " "),
          xlab = "median of out/in ratio", ylab = "decode_rate")
  }

  axis(1, at=1:bin, labels=round(tapply(df1.s$shift_ratio, df1.s$rank, median),2), las = 2)
  
  # axis(2, at=seq(0.6, 1, by=0.1), 
  #      labels=paste(seq(0.6, 1, by=0.1), sep=""))
  axis(2, at=seq(min(df1.s$avg_deco), max(df1.s$avg_deco), by=0.1),labels=paste(seq(min(df1.s$avg_deco), max(df1.s$avg_deco), by=0.1), sep=""))
  
  dev.off()
}







