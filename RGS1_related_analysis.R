rm(list=ls())

library(tidyverse)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(data.table)
library(dplyr)
library(ggpubr)

##############RGS1 expression across cancer types
col <- c("Normal"="#AC6C82", "Tumor"="#FFBC67", "Blood"="#DA727E")

a <- read.csv("./CD16low_RGS1.csv")
#a <- read.csv("./CD16high_RGS1.csv")

his_cell_num = a %>% filter((meta_tissue == "Tumor")|(meta_tissue == "Normal")|(meta_tissue == "Blood"))

his_cell_num = his_cell_num %>% separate(col = meta_histology, into = c("cancerType1", "cancerType2"), sep = "\\(")
his_cell_num <- his_cell_num %>% filter(!is.na(cancerType2))
his_cell_num = separate(data = his_cell_num, col = cancerType2, into = c("cancerType2", "cancerType3"), sep = "\\)")

cancernum = his_cell_num  %>% group_by(cancerType2,meta_tissue)

cancernum <- dplyr::summarise(cancernum,RGS1 = mean(RGS1))
cancernum<-data.frame(cancernum)


test <- dcast(cancernum,cancerType2~meta_tissue)
# test <- test[!is.na(test$Blood),]
# test[is.na(test)] <- 0
rownames(test) <- test$cancerType2
test <- test[,-1]

g1  <- pheatmap(
  test,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = "white",
  # cutree_cols = 2,
  # color = colorRampPalette(colors = c("white","red"))(30),
  cellwidth = 15, cellheight = 15
)
ggsave("CD16low_RGS1.pdf",g1, width = 4, height = 6)

g2  <- pheatmap(
  test,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = "white",
  # cutree_cols = 2,
  # color = colorRampPalette(colors = c("white","red"))(30),
  cellwidth = 15, cellheight = 15
)
ggsave("CD16high_RGS1.pdf",g2, width = 4, height = 6)

g <- ggarrange(g1,g2,heights=c(1, 1),ncol=2,common.legend = TRUE,legend="right",align = "v")

####RGS1 expression across tissues
a <- read.csv("CD16low_RGS1.csv")
#a <- read.csv("CD16high_RGS1.csv")

his_cell_num = a %>% filter((meta_tissue == "Tumor")|(meta_tissue == "Normal")|(meta_tissue == "Blood"))
his_cell_num = his_cell_num %>% separate(col = meta_histology, into = c("cancerType1", "cancerType2"), sep = "\\(")
his_cell_num <- his_cell_num %>% filter(!is.na(cancerType2))
his_cell_num = separate(data = his_cell_num, col = cancerType2, into = c("cancerType2", "cancerType3"), sep = "\\)")

cancernum = his_cell_num  %>% group_by(sampleID,meta_tissue)

cancernum <- dplyr::summarise(cancernum,RGS1 = mean(RGS1))
cancernum<-data.frame(cancernum)

col <- c("Normal"="#4974A4", "Tumor"="#B81316", "Blood"="#F29600")

sig_annotation = sig_test(cancernum,"meta_tissue","RGS1")
g2 <- cancernum %>% ggplot(aes(y=RGS1,x=meta_tissue,colour = meta_tissue))+
  # geom_violin( width=0.5,outlier.colour = NA, size=0.8)+
  geom_boxplot( width=0.5,outlier.colour = NA, size=0.8)+
  geom_jitter(size = 0.1,width = 0.12) + 
  ggsignif::geom_signif(data =sig_annotation %>% filter(p.value < 0.05),#  %>% filter(p.two.sided.fdr < 0.05),
                        mapping = aes(xmin = condition1, xmax = condition2, annotations = label, y_position = y),
                        manual = T,inherit.aes = F,textsize = 2.5)+
  theme_classic() +
  theme(axis.title.y=element_text(size=10),axis.text.y =element_text(size = 8),
        axis.text.x = element_text(angle = 60,size = 8 , hjust = 1, vjust = 1),
        strip.background = element_rect(colour = "white")) + 
  scale_fill_manual(values = col) +
  scale_colour_manual(values = col) +
  ylab("Means of log(RGS1)") +
  xlab("")


ggsave("CD16low_RGS1_boxplot_pancancer.pdf",g2, width = 3, height = 3)


######################RGS1_sensitivity_specificity

a <- read.csv("./fig2_bu_bili_0707.csv")
a
a$log10p <- -log10(a$cd16_p + 10^-301)
g <- ggplot(data = a,aes(x=cd16_score_ori,y=mean_cd16num_bi))+ 
  geom_point(aes( size = log10p),color='#CC99CC')+ #
  # scale_size_continuous(range = c(0,5))+
  # scale_color_gradient(low = "grey",high = "red") +
  # scale_color_manual(values = col) + 
  cowplot::theme_cowplot() +
  theme(axis.title.y=element_text(size=12),axis.text.y =element_text(size = 12),
        axis.text.x = element_text(size = 12),
        strip.background = element_rect(colour = "white")) + 
  ylab("Specifity score") +
  xlab("Sensitivity score") +
  geom_text_repel(size = 4,aes(label = X),color="black",max.overlaps=20)
g
ggsave("./permutation_prop_cd16_0707.pdf",g,width = 6, height = 5)

a$log10p <- -log10(a$cd56_p + 10^-301)
g <- ggplot(data = a,aes(x=cd56_score_ori,y=mean_cd56num_bi))+ 
  geom_point(aes( size = log10p),color='#CC99CC')+ #
  # scale_size_continuous(range = c(0,3))+
  # scale_color_gradient(low = "grey",high = "red") +
  # scale_color_manual(values = col) + 
  cowplot::theme_cowplot() +
  theme(axis.title.y=element_text(size=12),axis.text.y =element_text(size = 12),
        axis.text.x = element_text(size = 12),
        strip.background = element_rect(colour = "white")) + 
  ylab("Cell Type") +
  xlab("") +
  geom_text_repel(size = 4,aes(label = X),color="black",max.overlaps=20)
g
ggsave("./permutation_prop_cd56_0707.pdf",g,width = 6, height = 5)


######################
#heatmap
rm(list=ls())
CD16low <- read.csv("./CD16low_blood_tissue_DE_genes_0923.csv")
CD16low <- read.csv("./CD16high_blood_tissue_DE_genes_0923.csv")


CD16low = CD16low %>% filter((meta_tissue == "Tumor")|(meta_tissue == "Normal")|(meta_tissue == "Blood"))
# his_cell_num = DE %>% separate(col = meta_histology, into = c("cancerType1", "cancerType2"), sep = "\\(")
# his_cell_num[is.na(his_cell_num$cancerType2),"cancerType2"] = "Healthy donor"
# cancernum = separate(data = his_cell_num, col = cancerType2, into = c("cancerType2", "cancerType3"), sep = "\\)")
# 
# DE <- cancernum
# DE$Type <- "Patient"
# DE$Type[DE$cancerType2=="Healthy donor"] <- "Healthy donor"
# 
# DE <- DE[order(DE$Type,DE$cancerType2),]

Type <- CD16low$meta_tissue
# meta_histology <- DE$cancerType2
# majorType <- DE$majorType


CD16low <- CD16low[,c(4:(length(CD16low)))]

CD16low <- scale(CD16low,center=T,scale=T)
CD16low <- data.frame(CD16low)

CD16low[CD16low>2] = 2
CD16low[CD16low< -2] = -2
# CD16low$cancerType2 <-meta_histology
# CD16low$Type <- Type

col.heat <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))
cols <- col.heat(100)
# breaks <- seq(min(tmp), max(tmp), length = 100)
breaks <-seq(-2, 2, length = 100)
col1 = colorRamp2(breaks, cols)


exp_rng <- range(CD16low)
bks <- seq(exp_rng[1] - 0.1, exp_rng[2] + 0.1, by = 0.1)
col2 = colorRamp2(bks, blue2green2red(length(bks)))
  
pdf(paste0("CD16low_blood_tissue_DE_genes_0923.pdf"), width = 4, height = 4)
pdf(paste0("CD16high_blood_tissue_DE_genes_0923.pdf"), width = 4, height = 4)

ha <- rowAnnotation(
  Type = Type,
  # majorType=majorType,
  # Histology = meta_histology,
  col = list(Type = col <- c("Normal"="#4974A4", "Tumor"="#B81316", "Blood"="#F29600")
             # majorType = c("CD56highCD16low"="#023D93", "CD56lowCD16high"="#565656"),
             # Histology = col
  )
)

gene <- readRDS('gene_pos.rds')
gene_pos <- which(colnames(CD16low)%in%gene)
col_anno <- HeatmapAnnotation(gene=anno_mark(at=gene_pos,labels = gene),gp = gpar(fontsize=1))

g = Heatmap(as.matrix(CD16low), col = col1,
            show_row_names = FALSE,cluster_rows  = F,right_annotation = ha,
            show_column_names = FALSE,
            bottom_annotation =col_anno,
            # row_split = cluster_info,
            # column_names_gp = gpar( fontsize = 1),
            border = F)
print(g)
dev.off()

###############################################################################
###RGS1-other genes score

a <- read.csv("CD16low_all.csv", row.names = 1)
a <- read.csv("CD16high_all.csv",row.names = 1)

his_cell_num = a %>% filter((meta_tissue == "Tumor")|(meta_tissue == "Normal")|(meta_tissue == "Blood"))
# his_cell_num = his_cell_num %>% separate(col = meta_histology, into = c("cancerType1", "cancerType2"), sep = "\\(")
# his_cell_num <- his_cell_num %>% filter(!is.na(cancerType2))
# his_cell_num = separate(data = his_cell_num, col = cancerType2, into = c("cancerType2", "cancerType3"), sep = "\\)")

geneglist <- list("RGS1" = c("RGS1"),
                  "RGS1_CD69" = c("RGS1","CD69"),
                  "RGS1_ITGAE" = c("RGS1","ITGAE"),
                  "RGS1_CD69_ITGA1" = c("RGS1","CD69","ITGA1"),
                  "RGS1_CD69_ITGA1_ITGAE" = c("RGS1","CD69","ITGA1","ITGAE"),
                  "RGS1_CD69_ITGA1_ITGAE_ITGB1" = c("RGS1","CD69","ITGA1","ITGAE","ITGB1"),
                  "RGS1_CD69_ITGA1_ITGAE_ITGB1_CXCR6" = c("RGS1","CD69","ITGA1","ITGAE","ITGB1","CXCR6")
                  ) 


##AUC

meta <- his_cell_num[,c("meta_tissue","meta_histology","sampleID","batch","Majortype","celltype")]
his_cell_num <- t(his_cell_num[, !colnames(his_cell_num) %in% c("meta_tissue","meta_histology","sampleID","batch","Majortype","celltype")])
exprMatrix <- as.matrix(his_cell_num)
exprMatrix <- as(exprMatrix, "dgCMatrix")

cells_rankings <- AUCell_buildRankings(exprMatrix, plotStats=FALSE)
cells_AUC <- AUCell_calcAUC(geneglist, cells_rankings)

aucs <- getAUC(cells_AUC)

aucs_t <- data.frame(t(aucs))
aucs_t <- cbind(aucs_t,meta)

# write.csv(aucs_t,"./RGS1_other_AUC_CD16low.csv")
# write.csv(aucs_t,"./RGS1_other_AUC_CD16high.csv")
write.csv(aucs_t,"./RGS1_other_AUC_CD16low_0314.csv")
write.csv(aucs_t,"./RGS1_other_AUC_CD16high_0314.csv")
