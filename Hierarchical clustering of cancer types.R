library(tidyverse)
library(ggplot2)
library(ggtree)
library(treeio)
library(ggsci)
library(cowplot)

meta <- read.csv("meta_fix.csv",row.names = 1)

#Hierarchical clustering for cancer types
his_cell_num = meta %>% filter((meta_tissue == "Tumor"))

his_cell_num = his_cell_num %>% separate(col = meta_histology_new, into = c("cancerType1", "cancerType2"), sep = "\\(")
his_cell_num <- his_cell_num %>% filter(!is.na(cancerType2))
his_cell_num = separate(data = his_cell_num, col = cancerType2, into = c("cancerType2", "cancerType3"), sep = "\\)")


cancernum = his_cell_num%>%  group_by(cancerType2, celltype) 
test = his_cell_num %>% group_by(cancerType2) 


cancernum <- summarise(cancernum,count = n())
cancernum<-data.frame(cancernum)


test <- data.frame(summarise(test,
                             count = n()))

cancernum <- merge(test,cancernum,by="cancerType2")
cancernum['proportion'] <- cancernum$count.y / cancernum$count.x
head(cancernum)


cancernum2 <- cancernum[,c("cancerType2","celltype","proportion")]
cancernum2 <- dcast(cancernum2,cancerType2~celltype)
rownames(cancernum2) <- cancernum2$cancerType2
cancernum2 <- cancernum2[,-1]
cancernum2[is.na(cancernum2)] <- 0


tree = hclust(vegan::vegdist(cancernum2, method = 'bray'), 
              method = 'average')

p1 = ggtree(tree) + geom_tiplab() +xlim(NA,2)
p1

g2 <- cancernum2 %>% ggplot(aes(y=proportion,x=cancerType2,fill=celltype,colour = celltype))+
  geom_bar(stat = "identity",width = 0.5,alpha=0.65) +
  theme_classic() +
  
  theme(axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line = element_blank())+
  
  scale_x_discrete(limits = rev(c("OV","BRCA", "LC","CRC","HNSCC","ICC","HCC","PRAD","GC","RC","THCA","PACA","ESCA","BCC","UCEC","SCC","NPC","NB","FTC","MM","Leukemia", "MELA")),) +

  scale_fill_manual(values = c(CD16low_col, CD16high_col)) +
  scale_colour_manual(values = c(CD16low_col, CD16high_col)) +
  
  coord_flip()+
  ylab("Frequency") +
  xlab("") 
g2


g <- ggdraw()+
  draw_plot(p1, 0, 0.05, 0.7, 0.95)+
  draw_plot(g2, 0.2, 0, 0.8, 1)
g

ggsave("Hierarchical_clustering_for_Bright_Dim_NK_cells.pdf", g, height = 5, width =10)





#################
#Hierarchical clustering for AML/ALL/CLL
his_cell_num = meta %>% filter((meta_tissue == "Tumor")) %>% filter(meta_histology_new == "(Leukemia)")

his_cell_num = his_cell_num %>% separate(col = meta_histology, into = c("cancerType1", "cancerType2"), sep = "\\(")
his_cell_num <- his_cell_num %>% filter(!is.na(cancerType2))
his_cell_num = separate(data = his_cell_num, col = cancerType2, into = c("cancerType2", "cancerType3"), sep = "\\)")


cancernum = his_cell_num %>%filter((Majortype == "CD56lowCD16high") | (Majortype == "CD56highCD16high"))%>%  group_by(cancerType2, celltype) #cancerType2,
test = his_cell_num %>%  filter((Majortype == "CD56lowCD16high") | (Majortype == "CD56highCD16high"))%>%group_by(cancerType2) #cancerType2,

cancernum <- summarise(cancernum,count = n())
cancernum<-data.frame(cancernum)


test <- data.frame(summarise(test,
                             count = n()))

cancernum <- merge(test,cancernum,by="cancerType2")
cancernum['proportion'] <- cancernum$count.y / cancernum$count.x
head(cancernum)


cancernum2 <- cancernum[,c("cancerType2","celltype","proportion")]
cancernum2 <- dcast(cancernum2,cancerType2~celltype)
rownames(cancernum2) <- cancernum2$cancerType2
cancernum2 <- cancernum2[,-1]
cancernum2[is.na(cancernum2)] <- 0

tree = hclust(vegan::vegdist(cancernum2, method = 'bray'), 
              method = 'average')

p1 = ggtree(tree) + geom_tiplab() +xlim(NA,2)
p1

g2 <- cancernum%>% ggplot(aes(y=proportion,x=cancerType2,fill=celltype,colour = celltype))+
  geom_bar(stat = "identity",width = 0.5,alpha=0.65) +
  theme_classic() +
  
  theme(axis.title.y=element_text(size=12),axis.text.y =element_text(size = 12),
        axis.text.x = element_text(angle = 90,size = 12 , hjust = 1, vjust = 1),
        strip.background = element_rect(colour = "white")) +
  theme(axis.ticks.y = element_blank(),
        axis.line = element_blank())+
  
  scale_x_discrete(limits = c("AML","ALL","CLL"),) +

  scale_fill_manual(values = CD16high_col) +
  scale_colour_manual(values = CD16high_col) +
  
  coord_flip()+
  ylab("Frequency") +
  xlab("") 
g2
ggsave("Leukemia_NK_cells.pdf", g2, height = 3, width =7)

