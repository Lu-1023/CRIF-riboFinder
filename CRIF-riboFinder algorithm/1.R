library(qvalue)
library(readr)
library(ggrepel)
setwd("~/Data/FS_ss/FS_celegans/codon_repeat/FS")
samples<-system("ls ~/Data/ribo_seq/celegans/20210906/noscreen/coverage_stat/*_RPF_count_width0.txt",intern=T)
#RPF aim to add the RPKM and shift_ratio
RPF <- function(sam,f)  {
  #sam:the path of the RPF_count.txt, f:frame shifting sorts<1 or 2>
  x<-gsub(".*/(.*)_RPF_count_width0.txt","\\1",sam)
  d <- read.table(sam, sep="\t", header = T)
  d$ORF_length<-d$ORF_length-200
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
  #or
  FF<-read.table(paste0(x,"/",x,"frame",f,"_information.txt"),header = T,sep="\t")
  #A<-read.table(paste0(x,"/",x,"frame",f,"_array.txt"),header = T,sep="\t")
  #FF<-read.table(paste0(x,"/",x,"frame",f,"_hsc_information.txt"),header = T,sep="\t")
  total<-merge(FF,df1[,c(1,2,13)],by.x = "id",by.y="Gene",all.x = T)
  total$padjust<-p.adjust(total$p.value,method = "fdr",n=nrow(total))
  total_n<-total
  total<-na.omit(total)
  total$padjust<-p.adjust(total$p.value,method = "fdr",n=nrow(total))
  write.csv(total,paste0(x,"/",x,"frame",f,"_information_exp.csv"),row.names=FALSE)
  
  #A$id<-paste(A$id,A$start,sep="_")
  #A<-merge(total,A[,c(1,6,7)],by.x="id",by.y="id")
  #write.csv(A,paste0(x,"/",x,"frame",f,"_array.csv"),row.names=FALSE)
  if (f==1){
    total$color<-ifelse(total$padjust<0.05 & total$foldchange >1.5,"salmon","grey")
    color<-c("grey","salmon")
  }else{
    total$color<-ifelse(total$padjust<0.05 & total$foldchange >1.5,"steelblue3","grey")
    color<-c("grey","steelblue3")
  }
  #total<-total[total$foldchange>1.5,]
  exp_gene<-"" #c("EFNA3_NM_004952_529","KMT2D_NM_003482_11368","CPSF6_NM_001300947_1030")
  #total$id<-paste(total$id,total$start,sep="_")
  
  #png(filename = paste0(x,"/",x,"_fsgenes",f,"_volcano.png"),width = 1000,height =1000,bg="white")
  g<-ggplot(data=total,aes(x=log2(foldchange),y=-log10(padjust),color=color))+
    geom_point(size=3,shape=16)+
    labs(title = paste("algorism in ",f," frame"),y="-log10(p.adjust)",x= "log2(foldchange)")+
    scale_color_manual(values =color)+
    theme(plot.title = element_text(face = "bold",hjust=0.5,size = 40),
          axis.text.x = element_text(face="italic", color="grey20", size=30),
          axis.text.y = element_text( color="grey20", size=30),panel.background = element_rect(fill = "transparent",colour = "black",size=2.5),
          axis.title = element_text(size = 35,family = "Arial",color="grey20"),legend.text=element_text(size=35),
          legend.title=element_blank())+
    guides(fill=guide_legend(title=NULL))+
    scale_y_continuous(breaks = seq(0,25,3),limit=c(0,25))+
    scale_x_continuous(breaks = seq(-(round(max(log2(total$foldchange))+1)),round(max(log2(total$foldchange)))+1,2),
                       limit=c(-round(max(log2(total$foldchange))),round(max(log2(total$foldchange)))))+
    geom_hline(yintercept=-log10(0.05),lty=4,col="black",lwd=0.6)+
    geom_vline(xintercept=log2(1.5),lty=4,col="black",lwd=0.6)+
    geom_text_repel(data = subset(total,total$id %in% exp_gene),aes(label=id),size=10,box.padding = unit(10, "lines"),segment.color = "black")
  ggsave(filename = paste0(x,"/",x,"_fsgenes",f,"_volcano.eps"),plot=g,width = 1000,height =1000,bg="white",units = "px",dpi = 72)                  
  #print(g)
  #dev.off()
  #screen genes based on p-adjust value and foldchange
  total$color<-NULL
  total_fs<-total[total$padjust<0.05 &total$foldchange>1.5,]
  write.csv(total_fs,paste0(x,"/",x,"frame",f,"fs_fc1.csv"),row.names=FALSE)
  write.table(total_fs,paste0(x,"/",x,"frame",f,"fs_fc1.txt"),quote=FALSE,sep="\t",row.names=FALSE)
  total_fs$fs_sort<-rep(f,nrow(total_fs))
  write_csv(total_fs,paste0(x,"/",x,"fs_fc1.csv"),append = T,col_names = T)
  write.table(total[total$foldchange<1,],paste0(x,"/",x,"frame",f,"nc1.txt"),quote=FALSE,sep="\t",row.names=FALSE)
  total$fs_sort<-rep(f,nrow(total))
  write_csv(total,paste0(x,"/",x,"fs_expression.csv"),append = T,col_names = T)
}

for (i in samples) {
  RPF(i,1)
  RPF(i,2)
}
