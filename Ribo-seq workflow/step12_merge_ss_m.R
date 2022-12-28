#the script is used for get the start and stop counts in ribo-seq and compared with median

#1 Longest genes ,median and RPKM 
setwd("~/Data/ribo_seq/mouse_melanoma/1nonmerge/start_stop")
samples<-system("ls *.txt",intern = T)
filelist<-system("ls ~/Data/ribo_seq/hg38/20201207/coverage_stat/*RPF_count.txt",intern = T)

#RPF aim to add the RPKM and shift_ratio
RPF <- function(sam)  {
  d <- read.table(sam, sep="\t", header = T, stringsAsFactors = F)
  # calculate the number codon that are covered by at least 1 read 
  d$inframe_codon_cover_number = ceiling(d$inframe_cov/100 * (d$ORF_length - 200)/3 )
  # note  hist(d$inframe_codon_cover_number)
  d$RPF_total = apply(d[,c("inframe_count", "shift_plus_count", "shift_minus_count")], 1, sum)
  d$RPF_RPKM = d$RPF_total*1000000*1000/sum(d$RPF_total)/(d$ORF_length-200) 
  dfe<-d
  rm(d)
  dfe$shift_plus_ratio = dfe$shift_plus_count/dfe$RPF_total
  dfe$shift_minus_ratio = dfe$shift_minus_count/dfe$RPF_total
  dfe$shift_ratio = (dfe$shift_minus_count + dfe$shift_plus_count)/dfe$RPF_total
  df1 <- dfe[order(dfe$shift_ratio),]
  df1<-transform(df1,num_codon=(ORF_length-200)/3)
  df1<-transform(df1,cov=inframe_cov+shift_plus_cov+shift_minus_cov)
  #df1<-df1[df1$inframe_cov>5 & df1$RPF_RPKM>30 & df1$cov*df1$num_codon*0.01 > 20,] 
  df1<-df1[,c("Gene","RPF_RPKM","inframe_cov")]
  return(df1)
  
}

for ( x in seq(1,length(samples)) ) {
  ss<-read.table(samples[x],header = T,sep="\t")
  y<-gsub("(.*)_start_stop_counts_mRNA.txt","\\1",samples[x])
  M<-read.table(paste0("../decode/codon_median_",y,"_M_20_frame0.txt"),header = T,sep="\t")
  all<-merge(ss,M,by.x="gene_id",by.y = "X")mm10_fsYES_out_codon_same.txt
  all$median<-as.numeric(all$median)
  all<-transform(all,start2m=start_counts/median)
  all<-transform(all,stop2m=stop_counts/median)
  all<-all[order(all$stop_counts,decreasing = T),]
  rpf<-RPF(filelist[x])
  all<-merge(all,rpf,by.x="gene_id",by.y="Gene")
  write.csv(all,paste0(y,"start_stop_ratio_longest.csv"),row.names = F)
}


library(ggplot2)

#2 genome wide and median

setwd("~/Data/ribo_seq/chlamydomonas/start_stop")
samples<-system("ls *.txt",intern = T)
for ( x in samples ) {
  ss<-read.table(x,header = T,sep="\t")
  y<-gsub("(.*)_start_stop_counts_mRNA.txt","\\1",x)
  M<-read.table(paste0("../decode/codon_median_",y,"_M_20_frame0.txt"),header = T,sep="\t")
  M<-M[M$effective_codon!=0,]
  M$gene<-gsub("(.*)_.*","\\1",M$X)
  ss$gene<-gsub("(.*)_.*","\\1",ss$gene_id)
  all<-merge(ss,M,by.x="gene",by.y = "gene")
  all$median<-as.numeric(all$median)
  all<-transform(all,start2m=start_counts/median)
  all<-transform(all,stop2m=stop_counts/median)
  all<-all[order(all$stop_counts,decreasing = T),]
  all<-all[,-c(1,8)]
  write.csv(all,paste0(y,"start_stop_ratio_genomewide.csv"),row.names = F)
}
