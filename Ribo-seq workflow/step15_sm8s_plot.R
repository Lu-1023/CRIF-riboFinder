
get.max<-function(x){
  G<-read.csv(x,header = T)
  G<-transform(G,wb=gsub("(.*)_NM_.*","\\1",G$gene_id))
  Max<-aggregate(data=G,G$start_counts~wb,max)
  Max_G<-merge(G,Max,by.x="wb",by.y="wb")
  Max_G<-Max_G[Max_G$effective_codon*3/Max_G$total > 0.2 & Max_G$effective_codon>50,]
  Max_G<-unique(Max_G[,c(1,13,7)])
  return(Max_G)
}
setwd("~/Data/ribo_seq/mouse_melanoma/start_stop")
wide_list<-system("ls ~/Data/ribo_seq/mouse_melanoma/start_stop/*genomewide.csv",intern = T)

setwd("~/Data/ribo_seq/mouse_melanoma/sms_stat")
filelist<-list.files(path = "./",pattern = "RPF_distribution_meta.txt")
All<-data.frame()
#x is bin file ,z is the genomewide file ,nor is the normalize method (median or RNA RPKM)
get.bin<-function(x,z){
  a<-read.csv(x,header = F,sep="\t")
  #colnames(a)<-c("gene_ID","second","m1","m2","m3","m4","m5","m6","m7","m8","stop")
  a<-transform(a,V12=NULL)
  a<-transform(a,V43=NULL)
  a<-transform(a,wb=gsub("(.*)_NM_.*","\\1",a$V1))
  rMax_G<-get.max(z)
  n<-merge(rMax_G,a,by.x="wb",by.y="wb")
  n<-transform(n,V1=NULL)
  n_col<-colnames(n)
  n<-n[,c(n_col[1],n_col[4:13],n_col[2],n_col[14:43],n_col[3],n_col[44:53])]
  y<-gsub("(.*)_RPF.*","codon_median_\\1_M_20_frame0.txt",x)
  b<-read.table(paste0("../decode/",y),header=T,row.names = NULL)
  b<-b[,-c(3,4)]
  colnames(b)<-c("gene","median")
  b$median<-as.numeric(b$median)
  b<-transform(b,wb=gsub("(.*)_NM_.*","\\1",b$gene))
  b<-na.omit(b)
  
  #sometimes ,you can set the threshold to screen the data so that get a wonderful result
  #b<-b[b$median>1,]
  c<-merge(n,b,by.x="wb",by.y="wb")
  c<-transform(c,gene=NULL)
  d<-sweep(c[,2:53],1,c[,54],"/")
  rownames(d)<-c$gene_id
  e<-apply(d, 2, function(x) sum(na.omit(x))/length(na.omit(x)) )
  All<<-rbind(All,e)
}
for (i in seq(1,length(filelist))){
    get.bin(filelist[i],wide_list[i])
  
}
colnames(All)<-c(seq(0,51))
#All$ini<-c(0,0,0,0,0)
All<-t(All)
#colnames(All)<-c("recovering 15min","recovering 30min","recovering 60min","recovering 90min","non-treated")
colnames(All)<-c("wt-1","wt-2","wt-3","shELP3-1","shELP3-2","shELP3-3")
library(reshape2)
All<-transform(All,pos=c(seq(0,51)))
All<-melt(All,id.vars = "pos",value.name = "accumulation",variable.name = "sample")
png(filename = "Translation Initiation Pausing in 6 Samples.png",width =5000,height = 3000,bg="white")
ggplot(data=All,aes(x=pos,y=accumulation,group=sample,color=sample))+
  geom_line(size=5)+
  labs(title = "Translation Initiation Pausing in 6 Samples",y="PAUSING",x= "POSTION" )+
  scale_colour_manual(values =c("cornflowerblue","deeppink3","seagreen4","sienna1","wheat","yellow"),labels=c("wt-1","wt-2","wt-3","shELP3-1","shELP3-2","shELP3-3"),guide_legend(title="SAMPLES"))+
  theme(plot.title = element_text(face = "bold",hjust=0.5,size = 100),plot.subtitle = element_text(face = "bold",hjust=0.5,size =75),
        axis.text.x = element_text(color="grey20", size=70,angle=35),
        axis.text.y = element_text(face="italic", color="grey20", size=70),legend.title=element_text(size=90),legend.key.size = unit(60,"pt"),
        axis.title = element_text(size = 80,family = "Arial",color="grey20"),legend.text=element_text(size=90),
        panel.background = element_rect(fill = "transparent"),legend.position = "top")+
  scale_x_continuous(labels=NULL)+
  guides(fill=guide_legend(title="SAMPLES"))+
  annotate("rect",xmin = 0,xmax = 9,ymin = -0.35,ymax = -0.25,alpha=0.99,colour="black")+
  annotate("rect",xmin = 9,xmax = 41,ymin = -0.45,ymax = -0.18,alpha=0.2,colour="black")+
  annotate("rect",xmin = 41,xmax = 51,ymin = -0.35,ymax = -0.25,alpha=0.99,colour="black")+
  annotate("text", x = 10, y = -0.55, label = "2nd codon",size=20,color="orangered4")+
  annotate("text", x = 9, y = -0.85, label = "START\nCODON",size=25,color="grey47",fontface="bold")+
  annotate("text", x = 41, y = -0.85, label = "STOP\nCODON",size=25,color="grey47",fontface="bold")+
  geom_vline(xintercept = 10,linetype="dotted",size=3,color="orangered4")+
  geom_vline(xintercept = c(9,41),linetype="dotted",size=3,color="grey47")
  
dev.off()

#
setwd("~/Data/ribo_seq/mouse_melanoma/1nonmerge/sms_stat")
filelist<-list.files(path = "./",pattern = "RPF_distribution_meta.txt")
get.bin<-function(x){
  a<-read.csv(x,header = F,sep="\t")
  #colnames(a)<-c("gene_ID","second","m1","m2","m3","m4","m5","m6","m7","m8","stop")
  a<-transform(a,wb=gsub("(.*)_NM_.*","\\1",a$V1))
  y<-gsub("(.*)_RPF.*","codon_median_\\1_M_20_frame0.txt",x)
  b<-read.table(paste0("../decode/",y),header=T,row.names = NULL)
  b<-b[,-c(3,4)]
  colnames(b)<-c("gene","median")
  b$median<-as.numeric(b$median)
  b<-transform(b,wb=gsub("(.*)_NM_.*","\\1",b$gene))
  b<-na.omit(b)
  
  #sometimes ,you can set the threshold to screen the data so that get a wonderful result
  #b<-b[b$median>1,]
  c<-merge(a,b,by.x="wb",by.y="wb")
  c<-transform(c,gene=NULL)
  c<-transform(c,wb=NULL)
  d<-sweep(c[,2:53],1,c[,54],"/")
  rownames(d)<-c$V1
  e<-apply(d, 2, function(x) sum(na.omit(x))/length(na.omit(x)) )
  All<<-rbind(All,e)
}
All<-data.frame()
for (i in seq(1,length(filelist))){
  get.bin(filelist[i])
  
}
