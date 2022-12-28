#the script is used for dotploting to interpret the correlation between the expression and frameshift

setwd("~/Data/ribo_seq/chlamydomonas/coverage_stat")
library(ggplot2)
samples<-system("ls ~/Data/ribo_seq/chlamydomonas/coverage_stat/*RPF_count.txt",intern = T)

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
  #df1<-df1[,c("Gene","RPF_RPKM","shift_ratio")]
  return(df1)
 
}

changeFold<-function(x,y,name_x,name_y){
  a<-RPF(x)
  b<-RPF(y)
  a<-transform(a,RNA_seq_cov=NULL,RNA_seq_count=NULL)
  b<-transform(b,RNA_seq_cov=NULL,RNA_seq_count=NULL)
  all<-merge(a,b,by.x = "Gene",by.y = "Gene")
  all<-transform(all,RPKM_Fc=RPF_RPKM.y/RPF_RPKM.x,ratio_Fc=shift_ratio.y/shift_ratio.x)
  all<-na.omit(all)
  all<-all[all$ratio_Fc!=Inf,]
  all<-all[log2(all$RPKM_Fc)>0 | 0>log2(all$RPKM_Fc),]
  all<-all[all$ratio_Fc!=0,]
  #ggplot(data = all,aes(x=RPKM_Fc,y=ratio_Fc))+geom_point(size=3,shape=21)
  ggplot(data = all,aes(x=log2(RPKM_Fc),y=log2(ratio_Fc)))+geom_point(size=2,shape=16,color="darkolivegreen3",alpha=0.8)+
    labs(title = paste(name_x,"to",name_y),x="changes of gene expression(log2)",y="changes of frame-shift ratio(log2)")+
    theme(plot.title = element_text(face = "bold",hjust=0.5,size = 20),
          axis.text.x = element_text(face="italic", color="grey20", size=15),
          axis.text.y = element_text(size=15),
          axis.title = element_text(size = 15),
          panel.background = element_rect(fill = "transparent",colour = "gray65"),legend.position = c(0.9,0.9),)+
    annotate("text", x=4.5, y=3.5, label=c(paste("R =",round(cor(log2(all$RPKM_Fc),log2(all$ratio_Fc)),4))), family="serif", colour="lightsteelblue4", size=8)+
    annotate("text", x=4.5, y=3, label=c(paste("n =",nrow(all))), family="serif", colour="lightsteelblue4", size=8)
    
}


 # draw boxplot 
  jpeg(filename = paste(gsub("/home/song/Data/ribo_seq/chlamydomonas/coverage_stat/RiboSeq_(.*)_RPF_.*","\\1",sam), " dotplot.jpeg", sep=""),
       width = 350, height = 450, units = "px", quality = 100)
  
  par(mfrow=c(1,1))
 
  ggplot(data = df1,aes(x=log(RPF_RPKM),y=shift_ratio))+geom_point(size=3,shape=21)

  dev.off()
  
a<-RPF(samples[4])
b<-RPF(samples[5])
name_x<-gsub("/home/song/Data/ribo_seq/chlamydomonas/coverage_stat/RiboSeq_(.*)_RPF_.*","\\1",samples[4])
name_y<-gsub("/home/song/Data/ribo_seq/chlamydomonas/coverage_stat/RiboSeq_(.*)_RPF_.*","\\1",samples[5])
