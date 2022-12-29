##step1
library(ggplot2)
bg <- read.table("hom_dnds_bg.csv", sep = "\t")
cr = read.table("hom_dnds_cr.csv", sep = "\t")

crif <- read.csv("CRIF_total.csv")
crif$fid = paste(crif$id, crif$RNAID,sep = "_")

crif1 = merge(crif, cr, by.x = "fid", by.y = "V1")

bg$dns = (bg$V4+1)/(bg$V5+1)
cr$dns = (cr$V4+1)/(cr$V5+1)
crif1$dns = (crif1$V4+1)/(crif1$V5+1)

threshold = 98
d1 = crif1$dns[crif1$V3>threshold & crif1$CRIF_sort == "-1 CRIF"]
d1
write.table (d1, file ="hom_Riboseq_minus.csv", sep = ',')

boxplot(log2(bg$dns[bg$V3 >threshold ]), log2(cr$dns[cr$V3>threshold]), log2(crif1$dns[crif1$V3>threshold]),
        log2(crif1$dns[crif1$V3 > threshold & crif1$CRIF_sort == "+1 CRIF"]), log2(crif1$dns[crif1$V3>threshold & crif1$CRIF_sort == "-1 CRIF"]), notch = T )

##step2
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
piptide_file <- read.csv("hg38_mus_box.csv")
piptide_file
piptide_file <- data.frame(Type = piptide_file$Type, 
                           Peptide_len = as.numeric(piptide_file$Peptide_len), 
                           stringsAsFactors = FALSE)
names(piptide_file)<-c("Type","Peptide_len")
my_comparisons = list(c("All", "Homologs")
                      ,c("All", "C&H")
                      ,c("Homologs", "C&H")
)
cols <- c("black")
piptide_file$Type <- factor(piptide_file$Type,  levels = c("All", "Homologs","C&H"))
options(scipen = 200)
c <- ggplot(piptide_file, aes(x=Type, y=log10(Peptide_len), color = "black"))+
  geom_boxplot()+
  theme_bw()+
  theme(panel.grid = element_blank(), legend.position = c('none'), 
        axis.text.x = element_text(size = 20, face = "bold.italic"), 
        axis.text.y = element_text(size = 20, face = "bold.italic"), 
        axis.title.y = element_text(size = 20))+
  stat_compare_means(comparisons = my_comparisons, 
                     method = "t.test", 
                     #label= "p.signif",
  )+
  xlab("")
ggplot2::ggsave('hg38_mus_box.pdf', c, width = 10, 
                height = 15, limitsize = FALSE)