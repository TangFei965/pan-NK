###3. boxplot-violin
panC.freq.all.colSet.list <- readRDS("./panC.freq.all.colSet.list.rds")
col <- panC.freq.all.colSet.list[["cancerType"]]

#Fig1E
meta <- read.csv("./meta_fix.csv",row.names = 1)
col <- c("Normal"="#4974A4", "Tumor"="#B81316", "Blood"="#F29600")

his_cell_num = meta %>% filter((meta_tissue == "Tumor")|(meta_tissue == "Normal")|(meta_tissue == "Blood"))

his_cell_num = his_cell_num %>% separate(col = meta_histology_new, into = c("cancerType1", "cancerType2"), sep = "\\(")
his_cell_num <- his_cell_num %>% filter(!is.na(cancerType2))
his_cell_num = separate(data = his_cell_num, col = cancerType2, into = c("cancerType2", "cancerType3"), sep = "\\)")


cancernum = his_cell_num %>%filter((Majortype == "CD56lowCD16high") | (Majortype == "CD56highCD16high"))%>%  group_by(batch,meta_tissue, celltype) #cancerType2,
test = his_cell_num %>%  filter((Majortype == "CD56lowCD16high") | (Majortype == "CD56highCD16high"))%>%group_by(batch,meta_tissue) #cancerType2,

#Fig1F
cancernum = his_cell_num %>% filter(Majortype == "CD56highCD16low")%>%  group_by(batch,meta_tissue, celltype) #cancerType2,
test = his_cell_num %>%  filter(Majortype == "CD56highCD16low")%>%group_by(batch,meta_tissue) #cancerType2,


cancernum <- dplyr::summarise(cancernum,count = n())
cancernum<-data.frame(cancernum)


test <- data.frame(dplyr::summarise(test,
                                    count = n()))

cancernum <- merge(test,cancernum,by="batch")
cancernum <- cancernum[cancernum$meta_tissue.x == cancernum$meta_tissue.y,]
cancernum['proportion'] <- cancernum$count.y / cancernum$count.x
head(cancernum)

g2 <- cancernum %>% filter(count.x > 10) %>% ggplot(aes(y=proportion,x=celltype,colour = meta_tissue.x))+
  geom_violin( aes(fill=meta_tissue.x), width=0.8,outlier.colour = NA, size=0.5, scale="width")+
  geom_quasirandom(size = 0.05,width = 0.15) +
  stat_compare_means() +
  
  theme_classic() +
  theme(axis.title.y=element_text(size=12),axis.text.y =element_text(size = 12),
        axis.text.x = element_text(angle = 60,size = 8 , hjust = 1, vjust = 1),
        strip.background = element_rect(colour = "white"),
        axis.ticks =element_line(size=0.5),
        # aspect.ratio=1:2
  ) + 
  scale_fill_manual(values = col) +
  scale_colour_manual(values = col) +
  # ylab("Frequency of CD56dimCD16hi NK cells ") +
  ylab("Frequency of CD56highCD16low NK cells ") +
  xlab("") +
  facet_wrap(~meta_tissue.x, ncol = 3,scales = "free_y")
g2
ggsave("CD56highCD16low_proportion_violin_0311.pdf", g2, height = 4, width =7)
ggsave("CD56lowCD16high_proportion_violin_0311.pdf", g2, height = 3.5, width =12)


#Fig2B
his_cell_num = meta %>% filter((meta_tissue == "Tumor")|(meta_tissue == "Normal")) #|(meta_tissue == "Blood")

his_cell_num = his_cell_num %>% separate(col = meta_histology_new, into = c("cancerType1", "cancerType2"), sep = "\\(")
his_cell_num <- his_cell_num %>% filter(!is.na(cancerType2))
his_cell_num = separate(data = his_cell_num, col = cancerType2, into = c("cancerType2", "cancerType3"), sep = "\\)")

his_cell_num = his_cell_num[his_cell_num$cancerType2 %in% c("LC","RC","THCA","ESCA","BRCA","GC"),]

his_cell_num $celltype <- ifelse(as.vector(as.matrix(lapply(his_cell_num $celltype, function(x){grepl(x,pattern="CD16hi")}))),"CD16high","CD16low" )
cancernum = his_cell_num %>%  group_by(batch,meta_tissue,meta_histology, celltype) #cancerType2,
test = his_cell_num %>%group_by(batch,meta_tissue,meta_histology) #cancerType2,



cancernum <- dplyr::summarise(cancernum,count = n())
cancernum<-data.frame(cancernum)


test <- data.frame(dplyr::summarise(test,
                                    count = n()))

cancernum <- merge(test,cancernum,by="batch")
cancernum <- cancernum[cancernum$meta_tissue.x == cancernum$meta_tissue.y,]
cancernum['proportion'] <- cancernum$count.y / cancernum$count.x
head(cancernum)


col <- c("Normal"="#4974A4", "Tumor"="#B81316")
g <- 
  cancernum %>% ggplot(aes(y=proportion,x=celltype,colour = meta_tissue.x))+
  
  # geom_violin(aes(fill=meta_tissue.x), width=1,outlier.colour = NA, size=0.8,scale = "width")+
  geom_boxplot( fill="white", width=0.5,outlier.colour = NA, size=0.5, position = position_dodge(1))+
  geom_quasirandom(size = 0.5,width = 0.05, dodge.width = 1) +
  
  stat_compare_means(aes(group = meta_tissue.x), label = "p.format") +
  theme_classic() +
  
  theme(axis.title.y=element_text(size=12),
        axis.text.x = element_text(size = 12),
        legend.title = element_blank(),legend.position="none",
        strip.background = element_rect(colour = "white",fill = "white"),
        # legend.position="none",
        strip.text = element_text(size = 12, face = "bold" )) +  #,margins=TRUE)+ #axis.title.x =element_blank()
  scale_fill_manual(values = col) +
  scale_colour_manual(values = col) +
  ylab("Proportion in all NK cells") +
  xlab("") + facet_wrap(~meta_histology.x,nrow=2,scales = "free")
# print(g)
ggsave("proportion_box_majortype_NT_0311.pdf", g, height = 5, width =7)

#Fig3G
a <- read.csv("CD16low_RGS1.csv")
a <- read.csv("CD16high_RGS1.csv")

his_cell_num = a %>% filter((meta_tissue == "Tumor")|(meta_tissue == "Normal")|(meta_tissue == "Blood"))
his_cell_num = his_cell_num %>% separate(col = meta_histology, into = c("cancerType1", "cancerType2"), sep = "\\(")
his_cell_num <- his_cell_num %>% filter(!is.na(cancerType2))
his_cell_num = separate(data = his_cell_num, col = cancerType2, into = c("cancerType2", "cancerType3"), sep = "\\)")

cancernum = his_cell_num  %>% group_by(sampleID,meta_tissue)
#FigS4F
cancernum = his_cell_num  %>% filter(cancerType2=="CRC") %>% group_by(sampleID,meta_tissue)

cancernum <- dplyr::summarise(cancernum,RGS1 = mean(RGS1))
cancernum<-data.frame(cancernum)

col <- c("Normal"="#4974A4", "Tumor"="#B81316", "Blood"="#F29600")
sig_annotation = sig_test(cancernum,"meta_tissue","RGS1")
g1 <- cancernum %>% ggplot(aes(y=RGS1,x=meta_tissue,colour = meta_tissue))+
  
  geom_violin( aes(fill=meta_tissue), width=1,outlier.colour = NA, size=0.8, scale="width")+
  #geom_boxplot( fill="white", width=0.5,outlier.colour = NA, size=0.8)+
  geom_quasirandom(size = 0.2,width = 0.15) +
  
  # geom_jitter(size = 0.1,width = 0.12) + 
  ggsignif::geom_signif(data =sig_annotation %>% filter(p.value < 0.05),#  %>% filter(p.two.sided.fdr < 0.05),
                        mapping = aes(xmin = condition1, xmax = condition2, annotations = label, y_position = y),
                        manual = T,inherit.aes = F,textsize = 2.5)+
  theme_classic() +
  theme(axis.title.y=element_text(size=12),axis.text.y =element_text(size = 12),
        axis.text.x = element_text(size = 12 ),axis.ticks =element_line(size=0.5),
        strip.background = element_rect(colour = "white")) + 
  scale_fill_manual(values = col) +
  scale_colour_manual(values = col) +
  ylab("Means of log(RGS1)") +
  xlab("")+
  labs(title = "CD56dim CD16hi")
# labs(title = "CD56bright CD16lo")
# facet_wrap(~, ncol = 1,scales = "free_y")
g2

g <- g1 + g2
g

ggsave("RGS1_pancancer_0310.pdf", g, height = 3.5, width =7)
ggsave("RGS1_CRC_0131.pdf", g, height = 3.5, width =7)


#Fig5B
his_cell_num = meta %>% filter((meta_tissue == "Tumor"))

his_cell_num = his_cell_num %>% separate(col = meta_histology_new, into = c("cancerType1", "cancerType2"), sep = "\\(")
his_cell_num <- his_cell_num %>% filter(!is.na(cancerType2))
his_cell_num = separate(data = his_cell_num, col = cancerType2, into = c("cancerType2", "cancerType3"), sep = "\\)")

cancernum = his_cell_num %>% filter(meta_tissue=="Tumor") %>% group_by(batch,celltype,cancerType2)
cancernum <- dplyr::summarise(cancernum,count = n())
cancernum<-data.frame(cancernum)

test <- his_cell_num %>% filter(meta_tissue=="Tumor") %>%group_by(batch,cancerType2) 
test <- data.frame(dplyr::summarise(test,
                                    count = n()))

cancernum = his_cell_num %>%filter((Majortype == "CD56lowCD16high") | (Majortype == "CD56highCD16high"))%>%  filter(meta_tissue=="Tumor") %>% group_by(batch,celltype,cancerType2)#cancerType2,
test = his_cell_num %>%  filter((Majortype == "CD56lowCD16high")| (Majortype == "CD56highCD16high"))%>% filter(meta_tissue=="Tumor") %>%group_by(batch,cancerType2)  #cancerType2,

cancernum <- dplyr::summarise(cancernum,count = n())
cancernum<-data.frame(cancernum)
test <- data.frame(dplyr::summarise(test,
                                    count = n()))

cancernum <- merge(test,cancernum,by="batch")

cancernum <- cancernum[cancernum$cancerType2.x == cancernum$cancerType2.y,]
cancernum['proportion'] <- cancernum$count.y / cancernum$count.x
cancernum

cancernum <- cancernum %>% filter(celltype=="CD56dimCD16hi-c6-DNAJB1")
cancernum$cancerType2.x <- factor(cancernum$cancerType2.x,c("PRAD","GC","ESCA","HCC","RC","CRC","MELA","BRCA","PACA","NPC","LC","THCA","HNSCC","UCEC"))

g <- 
  cancernum %>% ggplot(aes(y=proportion,x=reorder(cancerType2.x, -proportion, FUN = median),colour = cancerType2.x))+
  geom_boxplot( fill="white", width=0.6,outlier.colour = NA, size=0.6, position = position_dodge(1))+
  geom_quasirandom(size = 0.5,width = 0.15) +
  
  theme_classic() +
  stat_compare_means() +
  theme(axis.title.y=element_text(size=12),
        axis.text.x = element_text(angle = 90,size = 12, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 12),
        legend.title = element_blank(),legend.position="none",
        strip.background = element_rect(colour = "white",fill = "white"),
        # legend.position="none",
        strip.text = element_text(size = 12, face = "bold" )) +
  scale_fill_manual(values = col) +
  scale_colour_manual(values = col) +
  ylab("Proportion of TaNK cells in CD56dimCD16hi NK cells") +
  xlab("")
g
ggsave("DNAJB1_prop.pdf", g, height = 2.5, width =4)

#FigS1D
meta <- read.csv("../NK_count.csv")
meta <- meta[!is.na(meta),]
his_cell_num = meta
his_cell_num = his_cell_num %>% separate(col = meta_histology, into = c("cancerType1", "cancerType2"), sep = "\\(")
cancernum <- his_cell_num %>% filter(!is.na(cancerType2))
cancernum = separate(data = cancernum, col = cancerType2, into = c("cancerType2", "cancerType3"), sep = "\\)")

cancernum <- cancernum[!is.na(cancernum$CD45),]

cancernum$proportion <- cancernum$cellID/cancernum$CD45


cancernum <- cancernum %>% filter(meta_tissue_new == "Tumor")
cancernum <- cancernum %>% filter(meta_tissue_new == "Normal")
cancernum <- cancernum %>% filter(meta_tissue_new == "Blood")


g <- ggplot(cancernum[cancernum$proportion<0.25,])+ 
  geom_violin(aes(x = reorder(cancerType2, proportion, FUN = median), y = proportion,colour = cancerType2),outlier.colour = NA,size=0.5,scale="width") + 
  #geom_boxplot(aes(x = reorder(cancerType2, proportion, FUN = median), y = proportion,colour = cancerType2),outlier.colour = NA,size=0.5,width=0.6) + 
  geom_quasirandom(aes(x = reorder(cancerType2, proportion, FUN = median), y = proportion,colour = cancerType2),size = 0.5,width = 0.15) +
  geom_hline(yintercept = median(cancernum$proportion[!is.na(cancernum$proportion)]),linetype = "dashed",colour="lightgrey",size=0.5) + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1,size = 12),
        axis.text.y = element_text(size = 12 ), axis.title.y = element_text(size = 12 ),
        panel.grid=element_line(colour=NA)) +#+
  scale_fill_manual(values = col) +
  
  # labs(title = "NK cells in tumor")+
  # labs(title = "NK cells in normal")+
  # labs(title = "NK cells in blood")+
  
  scale_colour_manual(values = col) +
  ylab("Propotioin in CD45+ cells") +
  xlab("")
print(g)

ggsave("NK_CD45_prop_tumor.pdf", g, height = 3, width =6)
ggsave("NK_CD45_prop_Normal.pdf", g, height = 3, width =4)
ggsave("NK_CD45_prop_Blood.pdf", g, height = 3, width =4)


#########################cellmarker

tmp <- read.csv("marker_gene_1026mean.csv",row.names = 1)
prop <- read.csv("marker_gene_1026prop.csv")
tmp <-scale(tmp)
tmp[tmp > 2] <- 2
tmp[tmp < -2] <- -2

tmp <- tmp[c('CD56brightCD16lo-c5-CREM',
             'CD56dimCD16hi-c4-NFKBIA',
             'CD56dimCD16hi-c6-DNAJB1',
             'CD56dimCD16hi-c7-NR4A3',
             
             'CD56dimCD16hi-c3-ZNF90',
             'CD56brightCD16hi',
             'CD56brightCD16lo-c2-IL7R-RGS1lo',
             'CD56brightCD16lo-c4-IL7R',
             
             'CD56brightCD16lo-c1-GZMH',
             'CD56dimCD16hi-c1-IL32',
             'CD56dimCD16hi-c2-CX3CR1',
             'CD56dimCD16hi-c5-MKI67',
             'CD56brightCD16lo-c3-CCL3',
             'CD56dimCD16hi-c8-KLRC2'),]

tmp = tmp[,rev(colnames(tmp))]

tmp <- melt(tmp)
prop <- melt(prop)


colnames(tmp) <- c("celltype","variable","expression")
tmp <- merge(tmp, prop,by=c("celltype","variable"))

tmp <- data.frame(tmp)
tmp$expression <- as.numeric(tmp$expression)


g <- ggplot(data = tmp)+ 
  geom_point(aes(x = celltype,y = variable, color = expression, size = value))+
  scale_size_area()+

  scale_color_gradientn(colors = c( '#1F78B4', '#FFFFB3', '#E31A1C')) +
  cowplot::theme_cowplot() +
  theme(axis.title.y=element_text(size=8),axis.text.y =element_text(size = 8),
        axis.text.x = element_text(size = 8,angle = 90),
        strip.background = element_rect(colour = "white")) + 
  ylab("Cell Type") +
  xlab("") + coord_flip()  
g
ggsave("marker_gene_1026.pdf", g, height = 4, width =14)