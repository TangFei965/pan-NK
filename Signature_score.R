### Calculation of signature score

#see details in Methods section for signature gene list
geneglist = readRDS('genelist.rds')

a <- read.csv("CD16high_all_log.csv",row.names = 1)
b <- read.csv("CD16low_all_log.csv",row.names = 1)

allcell = rbind(a,b)

meta <- allcell[,c("meta_tissue","meta_histology","sampleID","batch","Majortype","celltype")]
his_cell_num <- allcell[, !colnames(allcell) %in% c("meta_tissue","meta_histology","sampleID","batch","Majortype","celltype")]


meta = meta %>% separate(col = meta_histology, into = c("cancerType1", "cancerType2"), sep = "\\(")
meta[is.na(meta$cancerType2),"cancerType2"] <- "HD"
meta = separate(data = meta, col = cancerType2, into = c("cancerType2", "cancerType3"), sep = "\\)")



##AUC_score
exprMatrix <- as.matrix(his_cell_num)
exprMatrix <- as(exprMatrix, "dgCMatrix")

exprMatrix <-t(exprMatrix)

cells_rankings <- AUCell_buildRankings(exprMatrix, plotStats=FALSE)
cells_AUC <- AUCell_calcAUC(geneglist, cells_rankings)

aucs <- getAUC(cells_AUC)

aucs_t <- data.frame(t(aucs))
aucs_t$celltype <- meta[rownames(aucs_t), "celltype"]
aucs_t <- aucs_t %>% group_by(celltype) %>% dplyr::summarise_each(funs = mean)

write.csv(aucs_t ,"out/AUC.csv")


####plot
test1 <- aucs_t

test1$celltype <- factor(test1$celltype)

test1$Mc <- ifelse(grepl("CD56bright", test1$celltype),"CD56brightCD16lo","CD56dimCD16hi")
col = c('#003399', '#999999')
names(col) = c('CD56brightCD16lo', 'CD56dimCD16hi')



g1 <- test1  %>% ggplot() + geom_bar( aes(x =celltype ,y=inflammatory_score, fill=Mc),stat="identity",width = 0.8) + #,size = 3
  geom_hline(aes(yintercept=median(inflammatory_score)),linetype="dashed") + 
  ylab("Inflammatory score") +
  xlab("") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60,size = 0 , hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 12 ),
        panel.grid=element_line(colour=NA),
        panel.background = element_rect(fill = "#EDEDED"),
        axis.line = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = col) +
  scale_colour_manual(values = col) 
# coord_cartesian(ylim = c(0,25))
# facet_wrap(~Majortype,scales = "free",nrow=2)
g1

g2 <- test1  %>% ggplot() + geom_bar( aes(x =celltype ,y=Cytotoxicity, fill=Mc),stat="identity",width = 0.8) + #,size = 3
  geom_hline(aes(yintercept=median(Cytotoxicity)),linetype="dashed") + 
  ylab("Cytotoxicity score") +
  xlab("") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60,size = 0 , hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 12 ),
        panel.grid=element_line(colour=NA),
        panel.background = element_rect(fill = "#EDEDED"),
        axis.line = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = col) +
  scale_colour_manual(values = col) 

g3 <- test1  %>% ggplot() + geom_bar( aes(x =celltype ,y=stress_module, fill=Mc),stat="identity",width = 0.8) + #,size = 3
  geom_hline(aes(yintercept=median(stress_module)),linetype="dashed") + 
  ylab("Stress score") +
  xlab("") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60,size = 8 , hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 12 ),
        panel.grid=element_line(colour=NA),
        panel.background = element_rect(fill = "#EDEDED"),
        axis.line = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = col) +
  scale_colour_manual(values = col) 

# library(patchwork)
g <- g1 / g2 / g3
g
ggsave("./Fig1/Fig1_score_AUCell.pdf", g, height = 6, width =7)



############################################################################
#fig2_HLA_score
pdf("./Fig2/HLA_score_AUC.pdf", height = 6, width =7) 

cells_rankings <- AUCell_buildRankings(exprMatrix, plotStats=FALSE)
cells_AUC <- AUCell_calcAUC(geneglist, cells_rankings)

aucs <- getAUC(cells_AUC)

aucs_t <- data.frame(t(aucs))

test1 <- aucs_t


test1$celltype <- meta[rownames(test1),"celltype"]
test1$cancertype <- meta[rownames(test1),"cancerType2"]
test1$tissue <- meta[rownames(test1),"meta_tissue"]

test1 <- test1 %>% filter(tissue == "Tumor")

a = c("CD56brightCD16lo-c4-IL7R",
      "CD56brightCD16lo-c3-CCL3",  "CD56brightCD16lo-c5-CREM",  "CD56brightCD16hi",         
      "CD56dimCD16hi-c5-MKI67",    "CD56dimCD16hi-c4-NFKBIA",   "CD56dimCD16hi-c7-NR4A3",    
      "CD56dimCD16hi-c6-DNAJB1")

test1 = test1[test1$celltype %in% a,]


test1 <- test1 %>% group_by(cancertype, celltype) %>% dplyr::summarise_each(funs = mean)

for(j in colnames(test1)[c(8:11)]){
  print(j)
  tmp = test1[,c("cancertype" ,"celltype",j)]
  tmp = dcast(tmp,celltype~cancertype )
  
  rn = tmp$celltype
  tmp <- tmp[,-1]
  rownames(tmp) <- rn
  
  tmp <- tmp[,c("Leukemia","BRCA","CRC", "ESCA","GC",  "HCC", "HNSCC","LC", "MELA","NPC", "PACA","RC", "THCA")]
  tmp <-scale(tmp)
  tmp[tmp > 2] <- 2
  tmp[tmp < -2] <- -2
  # 
  # bk <- c(seq(0,0.5,by=0.01),seq(0.5001,1,by=0.01))
  g <- pheatmap(
    as.matrix(t(tmp)),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    gaps_col  = c(1,4),
    # breaks = bk,
    color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100), #length(bk)
    # color = colorRampPalette(c("navy", "white", "firebrick3"))(30),
    border_color = "white",
    # scale = "row",
    
    cellwidth = 10, cellheight = 10,main =paste0(i," ",j) , fontsize = 8
    
  )
  print(g)
}

dev.off()
