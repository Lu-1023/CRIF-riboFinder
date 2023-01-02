#plot line plot (10 codon upstream + 10 bin codon repeat box+ 10 codon downstream )
setwd("~/Data/FS_ss/FS_mouse/codon_repeat/FS_plot")
filelist<-system("ls -d */",intern = T)
mean_skip_zero<-function(x){
  a=0
  b=0
  for(i in x){
    if (i>0){
      a=a+1
      b=b+i
    } 
  }
  c=b/a
  return(b/a)
}
plot.bin.fs<-function(x,y,s){
  
  z<-gsub("(.*)/","\\1",x)
  
  if (s==1){
    fs<-read.table(paste0(x,z,"frame",y,"_1.txt"),header = F,sep="\t")
  }else if (s==2){
    fs<-read.table(paste0(x,z,"frame",y,"_2.txt"),header = F,sep="\t")
  }else{
    fs<-read.table(paste0(x,z,"frame",y,"_nc.txt"),header = F,sep="\t")
  }
  # 1 median normalized
  #media<-read.table(paste0(x,"codon_median_",z,"_M_0_frame",y,".txt"),header=T,sep="\t")
  #All<-merge(fs,media[,1:2],by.x="V1",by.y="gene")
  #All$median<-as.numeric(All$median)
  #All<-na.omit(All)
  
  #average normalized
  #fs$mean<-apply(fs[,6:35], 1, max)
  #fs$mean<-apply(fs[,6:35], 1, mean_skip_zero)
  fs$mean<-apply(fs[,6:35], 1, mean)
  ALL<-fs
  
  
  #merge script
  #ALL<-sweep(ALL[,6:45],1,ALL[,46],'/')
  ALL<-sweep(ALL[,6:35],1,ALL[,36],'/')
  ALL[is.na(ALL)]<-0
  #e<-apply(ALL[1:40], 2, function(x) sum(na.omit(x))/length(na.omit(x)) )
  e<-apply(ALL[1:30], 2, function(x) sum(na.omit(x))/length(na.omit(x)) )
  total<<-rbind(total,e)
  return(ALL)
}

library(ggalt)
frame<-2
co<-c("salmon","steelblue3")#"salmon" #
c<-""    #"median"
for (i in seq(1,length(filelist))){
  total<-data.frame()
  a<-plot.bin.fs(filelist[i],frame,1)
  b<-plot.bin.fs(filelist[i],frame,2)
  c<-plot.bin.fs(filelist[i],frame,"nc")
  #colnames(total)<-c(seq(-15,24))
  colnames(total)<-c(seq(-10,19))
  total<-t(total)
  total<-data.frame(total)
  #
  #colnames(total)<-c(seq(1,nrow(a)),seq(nrow(a)+1,nrow(b)+nrow(a)))
  #
  colnames(total)<-c("FS_R","FS","NC")
  library(reshape2)
  total$pos<-as.numeric(rownames(total))
  total<-melt(total,id.vars = "pos",value.name = "counts",variable.name = "sample")
  #
  #total$sample<-c(rep("FS",nrow(a)),rep("NC",nrow(b)))
  #png(filename = paste0(gsub("(.*)/","\\1",filelist[i]),"_the reads count for repeat box",frame,"_",c,".png"),width =5000,height = 2400,bg="white")
  g<-ggplot(data=total,aes(x=pos,y=counts,group=sample,color=sample))+
    geom_xspline(size=15)+
    labs(title = paste0("the Reads Counts for Repeat-box (Frame ",frame,") in ",filelist[i]),y="read counts\n(sum counts normalized)",x="")+
    scale_colour_manual(values =c(co,"grey73","seagreen4","sienna1","wheat"),labels=c("CRAF_r events","CRAF events","negative control"),guide_legend(title=""))+
    theme(plot.title = element_text(face = "bold",hjust=0.5,size = 100),plot.subtitle = element_text(face = "bold",hjust=0.5,size =75,family = "Arial"),
          axis.text.x = element_text(color="grey20", size=70,angle=35,vjust=10,family = "Arial"),
          axis.text.y = element_text(color="grey20", size=70,family = "Arial"),legend.title=element_text(size=90),legend.key.size = unit(60,"pt"),
          axis.title = element_text(size = 80,family = "Arial",color="grey20"),legend.text=element_text(size=90),
          panel.background = element_rect(fill = "transparent"),legend.position = "right",
          axis.line = element_line(colour = "black",size=5))+
    scale_x_continuous(breaks=seq(-10,20,10),expand = c(0,0),labels = NULL)+
    scale_y_continuous(expand = c(0,0))+
    guides(fill=guide_legend(title=""))+
    annotate("text", x = 0, y = min(total$counts)+0.007, label = "codon repeat start",size=35,color="black",family = "Arial")+
    annotate("text", x = 9, y = min(total$counts)+0.007, label = "hidden stop codon",size=35,color="black",family = "Arial")+
    geom_vline(xintercept = 0,linetype="dotted",size=5,color="grey47")+
    geom_vline(xintercept = 9,linetype="dotted",size=5,color="grey47")
  ggsave(filename = paste0(gsub("(.*)/","\\1",filelist[i]),"_the reads count for repeat box",frame,".eps"),plot=g,width =5000,height = 2400,bg="white",units = "px",dpi = 72,limitsize = FALSE)
  #print(g)
  #dev.off()
}
