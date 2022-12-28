#The script can get the translation efficience (FPKM.Ribo_seq/FPKM.RNA_seq)

#1 read RNA data
setwd("~/Data/ribo_seq/celegans/ballgown2")
filelist<-list.files(path = "./",pattern = "m")
filelist<-c(filelist[9:10],filelist[1:8])
read.rna<-function(x){
  phenotype<-data.frame(sample=x,treat="celegans")
  b<-ballgown(dataDir = "./",samples= x,pData=phenotype)
  g<-gexpr(b)
  rm(b)
  g<-transform(g,gene=rownames(g))
  return(g)
}
file_rpkm<-system("ls ~/Data/ribo_seq/celegans/start_stop/start_stop/*longest.csv",intern = T)
te.cal<-function(x,y){
  ribo<-read.csv(x,header = T)
  ribo<-ribo[,c(1,13)]
  ribo<-transform(ribo,gene=gsub("(.*)_.*","\\1",gene_id))
  total<-merge(y,ribo,by.x="gene",by.y="gene")
  total<-transform(total,te=RPF_RPKM/avg)
  total<-total[total$avg!=0,]
  sam<-gsub(".*stop/(.-.*)start.*csv","\\1",x)
  write.csv(total,paste0(sam,"_translation_efficience.csv"),row.names  = F)
}
cutoff<-"1to2" # you can set "1to1" or "1to2" , that means you have one Ribo_RPKM,but two RNA_RPKM
for (i in seq(1,length(file_rpkm))){
  if (cutoff=="1to2"){
    a1<-read.rna(filelist[2*i-1])
    a2<-read.rna(filelist[2*i])
    All<-merge(a1,a2,by.x="gene",by.y="gene")
    All$avg<-apply(All[,2:3],1,mean)
    te.cal(file_rpkm[i],All)
  }else{
    te.cal(file_rpkm[i],filelist[i])
  }
}

#correlation between pausing and translation efficiency
sta<-read.csv("2-sta_translation_efficience.csv",header = T)
r10<-read.csv("3-sr10_translation_efficience.csv",header = T)
TE<-merge(sta[,c(1,7)],r10[,c(1,7)],by.x="gene",by.y="gene")
TE<-transform(TE,fc_te=log2(te.y/te.x))
te_pausing<-merge(TE[,c(1,4)],point2[,c(1,10)],by.x="gene",by.y="gene")
ggplot(data = te_pausing,aes(x=pausing_c,y=fc_te))+geom_point(size=3,shape=16,color="darkolivegreen3",alpha=0.8)+
  labs(title = paste("starvation","to","recovering 10min"),x="Pausing Foldchange(log2)",y="Translation Efficiency(log2)")+
  theme(plot.title = element_text(face = "bold",hjust=0.5,size = 30),
        axis.text.x = element_text(face="italic", color="grey20", size=25),
        axis.text.y = element_text(size=25),
        axis.title = element_text(size = 25),
        panel.background = element_rect(fill = "transparent",colour = "gray65"),legend.position = c(0.9,0.9),)+
  annotate("text", x=-1.2, y=3.5, label=c(paste("R =",round(cor(te_pausing$fc_te,te_pausing$pausing_c),4))), family="serif", colour="lightsteelblue4", size=12)+
  annotate("text", x=-1.2, y=3, label=c(paste("n =",nrow(te_pausing))), family="serif", colour="lightsteelblue4", size=12)
