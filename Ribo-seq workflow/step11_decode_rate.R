#the script is used for calculating decode rate by RNA_Seq or median

setwd("~/Data/ribo_seq/celegans/decode")

#1 using RNA_seq(note:line20,)
detail<-read.table("detail.txt",header = F)
samples<-detail[,1]
for ( x in samples )    {
  
  number<-read.table(paste("codon_number_",x,"_M_20_frame0.txt",sep = ""),header=T,sep = "\t")
  reads<-read.table(paste("codon_reads_",x,"_M_20_frame0.txt",sep = ""),header=T,sep = "\t")
  name<-reads[,1]
  reads<-reads[,-c(1,2)]
  number<-number[,-1]
  number<-as.matrix(number)
  reads<-as.matrix(reads)
  #sum(reads) assume frame 0 total reads
  reads<-reads*1000000/sum(reads)
  n_value<-reads/number
  #Note!!! gene name you can use the (1)_(2)  1 or 2 ,better for symbol name . 
  rownames(n_value) = gsub("(.*)_.*","\\1",name)
  n_value<-transform(n_value,gene=gsub("(.*)_.*","\\1",name))
  
  phenotype<-data.frame(id=x,treat="treat")
  #Note!!! the ballgown file name
  b<-ballgown(dataDir = paste("~/Data/ribo_seq/chlamydomonas/RNA_seq/ballgown/",gsub(".*_(.*)","\\1",x),sep=""),samplePattern = "RNASeq",pData =phenotype)
  t<-texpr(b,meas = "all")
  rm(b)
  #using Ribo_seq FPKM to screen
  # phenotype1<-data.frame(id=y,treat="treat")
  # b1<-ballgown(dataDir = paste("~/Data/ribo_seq/chlamydomonas/fq/ballgown/",gsub(".*_(.*)","\\1",y),sep=""),samplePattern = "RiboSeq",pData =phenotype1)
  # t1<-texpr(b1,meas = "all")
  # rm(b1)
  # t1<-t1[,c(6,12)]
  # t1<-t1[t1[,2]>0,]
  #total<-merge(n_value,t1,by.x = "gene",by.y = "t_name")
  #procesing RNA-seq t
  t<-t[,c(6,12)]
  t<-t[t[,2]>1,]
  total<-merge(n_value,t,by.x = "gene",by.y = "t_name")
  #total<-merge(total,t,by.x = "gene",by.y = "t_name")
  #total<-total[total[,67]>0,]
  #resort the colomn of total
  # cols<-colnames(total)
  # new_cols<-c(cols[1:65],cols[67],cols[66])
  # total<-total[,new_cols]
  #calculate decode rate,or use command : dtotal<-total[,2:63]/total[,64]
  dtotal<-sweep(total[,2:65],1,total[,66],'/')
  rownames(dtotal)<-total$gene
  dtotal<-as.matrix(dtotal)
  
  # using apply(dtotal, 2, function(x) sum(na.omit(x))/length(na.omit(x)) )
  a<-0
  for (i in seq(1,64)) {
    avg<-sum(na.omit(dtotal[,i]))/length(na.omit(dtotal[,i]))
    print(avg)
    a<-append(a,avg)
  }
  a<-a[-1]
  dtotal<-rbind(dtotal,a)
  write.csv(dtotal,paste("RNASEQ_",x,".csv",sep=""))
}

#2 using median(note:line78,85,91)
detail<-read.table("detail.txt",header = F)
samples<-detail[,1]
for ( x in samples )    {
  
  number<-read.table(paste("codon_number_",x,"_M_20_frame0.txt",sep = ""),header=T,sep = "\t")
  reads<-read.table(paste("codon_reads_",x,"_M_20_frame0.txt",sep = ""),header=T,sep = "\t")
  m<-read.table(paste("codon_median_",x,"_M_20_frame0.txt",sep = ""),header=T,sep = "\t")
  name<-reads[,1]
  reads<-reads[,-c(1,2)]
  number<-number[,-1]
  m<-m[m[,3]!=0,]
  m<-transform(m,X=gsub("(.*)_(.*)","\\2",X))
  m<-transform(m,ratio=effective_codon/total_codon)
  # screening m using effecitve_codon or total_codon or ratio to get about 200 genes,here we select genes' effective_codon>50
  e_m<-m[m$ratio>0.1 & m$effective_codon>40,]
  number<-as.matrix(number)
  reads<-as.matrix(reads)
  reads<-reads*1000000/sum(reads)
  n_value<-reads/number
  
  rownames(n_value) = gsub("(.*)_(.*)","\\2",name)
  n_value<-transform(n_value,gene=rownames(n_value))
  e_m<-as.matrix(e_m)
  n_value<-merge(n_value,e_m,by.x="gene",by.y="X")
  etotal<-sweep(n_value[,2:65],1,as.numeric(n_value[,66]),'/')
  #note here!!! or 
  letotal<-etotal
  #letotal<-log2(etotal) 
  row.names(letotal)<-n_value[,1]
  letotal<-rbind(letotal,apply(letotal, 2, function(x) sum(na.omit(x))/length(na.omit(x)) ))
  write.csv(letotal,paste0("MEDIAN_",x,".csv"))
  
}

#3 calculate the genes number and genes_reads in RNASEQ_x.csv or MEDIAN_x.csv
detail<-read.table("detail.txt",header = F)
samples<-detail[,1]
for (x in samples) {
  #y<-gsub(".*_(.*)","RiboSeq_\\1",x)
  reads<-read.table(paste("codon_reads_",x,"_M_20_frame0.txt",sep = ""),header=T,sep = "\t")
  name<-gsub("(.*)_(.*)","\\2",reads[,1])
  reads<-reads[,-c(1,2)]
  reads<-as.matrix(reads)
  #sum(reads) assume frame 0 total reads
  reads<-reads*1000000/sum(reads)
  tem<-read.csv(paste("MEDIAN_",x,".csv",sep=""),header = T)
  tem<-tem[-nrow(tem),]
  reads<-transform(reads,gene=name)
  rownames(reads)<-name 
  gene_reads<-reads[tem[,1],]
  s<-sum(as.matrix(gene_reads[,1:64]))
  print(paste(x,nrow(tem),s,s/1000000))
  
}  

#4  merge the information of AA codon CAI decode_rate to a table
#first make i<-detail[1,1],run the loop ,second rownames(AVG) contains the 64 codons
detail<-read.table("detail.txt",header = F)
samples<-detail[-1,1]
all<-data.frame(row.names = rownames(AVG))
for (i in samples)   {
  #tem<-read.csv(paste("genelist/",i,".csv",sep=""),header = T,row.names = 1)
  tem<-read.csv(paste("MEDIAN_",i,".csv",sep=""),header = T,row.names = 1)
  AVG<-tem[nrow(tem),]
  AVG<-t(AVG)
  colnames(AVG)<-i  
  all<-cbind(all,AVG)
  
}

#5 translate code to AA
#BiocManager::install("seqinr")
#library(seqinr)
library(Biostrings)
AA<-""
for (i in seq(1,length(rownames(all)))) {
  code<-DNAString(rownames(all)[i])
  aa<-translate(code,no.init.codon = TRUE)
  AA<-append(AA,as.character(aa))
  
}
AA<-AA[-1]
all<-transform(all,AA=AA)
all<-transform(all,codon=rownames(all))
all<-all[order(all[,(ncol(all)-1)],all[,ncol(all)]),]
CAI<-read.table("CBI_ref.txt",header = T)
ALL<-merge(all,CAI,by.x="codon",by.y="codon")
ALL<-ALL[order(ALL["AA.x"],ALL["codon"]),]
ALL<-ALL[,-c(ncol(ALL),(ncol(ALL)-2))]
#resorts the colomns
cols<-colnames(ALL)
new_cols<-c(cols[1],cols[(ncol(ALL)-1):(ncol(ALL))],cols[2:(ncol(ALL)-2)])
ALL<-ALL[,new_cols]
write.csv(ALL,"total_DR.csv",row.names = F)
rownames(ALL)<-paste(ALL$codon,ALL$AA.x,sep = "_")
library(pheatmap)

#  eps quality is higher than jpeg jpeg(filename = width = 800,height = 1000,units = "px", bg="white",quality = 100)
jpeg(filename = "decoedrate.jpeg",width = 800,height = 1000,units = "px", bg="white",quality = 100)
#log2
h<-pheatmap(log2(as.matrix(ALL[4:64,4:ncol(ALL)])),cluster_cols = F,cluster_rows = F,cellwidth=20,cellheight=10,show_rownames = T,
            border_color=NA,color=colorRampPalette(c("dodgerblue","white","brown1"))(10),main="DECODE_RATE (log2))")
#plot directly
h<-pheatmap(as.matrix(ALL[4:64,4:ncol(ALL)]),cluster_cols = F,cluster_rows = F,cellwidth=20,cellheight=10,show_rownames = T,
            border_color=NA,color=colorRampPalette(c("dodgerblue","white","brown1"))(10),main="DECODE_RATE in celegans")

#title(main = "DECODE_RATE in Chlamydomonas",cex.main=5,font=3)
dev.off()
detail<-read.table("detail.txt",header = F)
samples<-detail[,1]
for (i in samples)   {
  cbiref(i)
}
cbiref = function(x) {
  
  df = read.csv(paste("MEDIAN_",x,".csv",sep=""), header = T,row.names = 1)
  decode_rate = df[nrow(df),]
  decode_rate<-decode_rate[order(decode_rate)]
  decode_rate<-t(decode_rate)
  ref<-data.frame(decode_rate=decode_rate)
  ref$rand<-c(rep(3,22),rep(2,21),rep(1,21))
  
  names(ref) = c("decode_rate", "rand")
  
  write.table(ref, paste0(x, "_DECODE_ref.txt"), sep = "\t", row.names = T, quote = F)
  
}


