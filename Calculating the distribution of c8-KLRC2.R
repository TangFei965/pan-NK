#Ro/e

meta <- read.csv("./meta_fix.csv",row.names = 1)

his_cell_num = meta %>% filter(meta_tissue == "Blood") 

his_cell_num = his_cell_num %>% separate(col = meta_histology_new, into = c("cancerType1", "cancerType2"), sep = "\\(")
his_cell_num <- his_cell_num %>% filter(!is.na(cancerType2))
his_cell_num = separate(data = his_cell_num, col = cancerType2, into = c("cancerType2", "cancerType3"), sep = "\\)")

cancernum = his_cell_num  %>% group_by(cancerType2,celltype)
cancernum <- dplyr::summarise(cancernum,count = n())
cancernum<-data.frame(cancernum)


test <- dcast(cancernum,cancerType2~celltype)
test[is.na(test)] <- 0
rownames(test)<-test$cancerType2
test<-test[,-1]
test <- as.matrix(test)

a =margin.table(test, 1)
b=margin.table(test, 2)
c=margin.table(test)
test = test/(outer(a,b,"*")/c)

tmp <- melt(test)
tmp <- tmp[order(tmp$Var1,tmp$value,tmp$Var2),]

tmp$Var1 <- factor(tmp$Var1,levels=as.character(unique(tmp$Var1)))
tmp$state <- ifelse(tmp$value>1,"Enrichment","Depletion")

# col = c("#AF2387","#003399")
# names(col) = c("Enrichment","Depletion")

g <- ggplot(data = tmp)+ 
  geom_point(aes(x = Var1,y = Var2, color = state, size = value))+
  scale_size_area()+
  # scale_color_gradient(low = "grey",high = "red") + 
  scale_color_manual(values = col) + 
  cowplot::theme_cowplot() +
  theme(axis.title.y=element_text(size=12),axis.text.y =element_text(size = 12),
        axis.text.x = element_text(size = 12),
        strip.background = element_rect(colour = "white")) + 
  ylab("Cell Type") +
  xlab("") 
g

tmp <- tmp[tmp$Var2=="CD56dimCD16hi-c8-KLRC2",]
colnames(tmp) <- c("cancerType2","celltype","RO/E","state")


g<-ggdotchart(tmp, x = "cancerType2", y = "RO/E",
              color = "cancerType2",                              # Color by groups
              palette =  col, # Custom color palette
              sorting = "descending",                       # Sort value in descending order
              add = "segments",                             # Add segments from y = 0 to dots
              add.params = list(color = "lightgray", size = 2), # Change segment color and size
              dot.size = 8,                                 # Large dot size
              # font.label = list(color = "white", size = 3,
              #                   vjust = 0.5),               # Adjust label parameters
              # rotate = TRUE,
              ggtheme = theme_pubr()                        # ggplot2 theme
)+ geom_hline(yintercept = 1, linetype = 2, color = "gray") #+ geom_hline(yintercept = -1, linetype = 2, color = "gray") 
g


ggsave(paste0("./figS7/","c8-KLRC2-ROE.pdf"),g, width = 4, height = 5)
