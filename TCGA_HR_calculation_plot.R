library(dplyr)
library(ggplot2)
library(pheatmap)
library(data.table)
library(circlize)
library(grDevices)
library(survminer)
library(survival)
library(RegParallel)
library(ggpubr)
library(corrplot)
library(patchwork)

rm(list = ls())


adf = read.table('TCGA_aucs.csv', header=T, stringsAsFactors=F, check.names=F, sep="\t")

dat = read.table("TCGA_surdf.csv", header=T, stringsAsFactors=F, check.names=F, sep=",")
dat = dat[,c("sample","X_PATIENT","cancer.type.abbreviation","age_at_initial_pathologic_diagnosis",
             "gender","race","ajcc_pathologic_tumor_stage",
             "OS","OS.time","sampleType")]
rownames(dat)=dat$sample
rownames(adf)=adf$sample
use=intersect(rownames(dat),rownames(adf))
dat = dat[use,]
adf = adf[use,]


dat = cbind(dat,adf)
dat = dat[!is.na(dat$OS.time),]
dat = dat[!is.na(dat$age_at_initial_pathologic_diagnosis),]
dat$age = dat$age_at_initial_pathologic_diagnosis
dat$cancerType = dat$cancer.type.abbreviation
dat = dat[!is.na(dat$ajcc_pathologic_tumor_stage),]
dat$stage = dat$ajcc_pathologic_tumor_stage
unique(dat$stage)
dat = dat[grepl("Stage I",dat$stage),]
unique(dat$stage)
dat$stage = gsub("[ABCD]$","",dat$stage)
unique(dat$stage)

### patch: merge COAD and READ as CRC; remove typing biased cancers
dat$cancerType = ifelse(dat$cancerType=="COAD"|dat$cancerType=="READ", "CRC", dat$cancerType)


dat = dat[dat$cancerType %in% c("BLCA","HNSC","PAAD","CRC","LIHC","BRCA","LUAD","ESCA","STAD",
                                'SKCM','UVM','KIRC','MESO'),]

dat$OS.time <- dat$OS.time / 365

dat$OS[dat$OS.time>7] <- 0
dat$OS.time[dat$OS.time>7] <- 7

dat$OS.use=dat$OS
dat$OS.time.use=dat$OS.time

clus=c("aucs")
dir="0.out"

do.cox = function(df, type, cancer){
  #print(unique(df$group))
  
  df$group = factor(df$group,levels = c('low','high'))
  cox=c()
  if (type=="reg"){
    if(length(levels(df$gender)>1)){
      cox = coxph(Surv(OS.time.use,OS.use) ~ group + gender + age + stage, data=df)
    }else{
      cox = coxph(Surv(OS.time.use,OS.use) ~ group + age + stage, data=df)
    }
  }else{
    cox = coxph(Surv(OS.time.use,OS.use) ~ group, data=df)
  }
  cox.summary = summary(cox)
  #print(cox)
  #print(cox.summary)
  nvar = length(unique(df$group)) - 1
  #
  HR = round(cox.summary$conf.int[1:nvar,1], 2)
  HR.lower = round(cox.summary$conf.int[1:nvar,3], 2)
  HR.upper = round(cox.summary$conf.int[1:nvar,4], 2)
  HR.range = sprintf("(%.1f-%.1f)", HR.lower, HR.upper)
  coef = cox.summary$coefficients[1:nvar,1]
  coef.se = cox.summary$coefficients[1:nvar,3]
  Pval = round(cox.summary$coefficients[1:nvar,5], 4)
  group = gsub("group","",rownames(cox.summary$conf.int)[1:nvar])
  return(data.frame(cancerType = cancer, group=group, HR=HR, HR.range=HR.range, coef=coef, coef.se=coef.se, Pval=Pval))
}

for (i in clus) {
  dir.create(sprintf("%s/HR1/%s",dir,i), F, T)
  
  types = c("reg",'noreg')
  for (type in types){
    cdat = data.frame(cancerType = 'new', group='new', HR=1, HR.range='(0-0)', coef=1, coef.se=1, Pval=1)
    for (cancer in c("panCan", unique(dat$cancerType))){
      pdat = dat
      if(cancer!="panCan"){
        pdat = dat[dat$cancerType==cancer,]
      }
      res.cut <- surv_cutpoint(pdat, time = "OS.time.use", 
                               event = "OS.use", 
                               variables = i, 
                               )
      res.cat <- surv_categorize(res.cut)
      pdat$group <- res.cat[,i] 
      pdat$group = factor(pdat$group)
      
      test = do.cox(pdat,type,cancer)
      cdat = rbind(cdat,test)
    }
    
    cdat = cdat[cdat$cancerType!='new',]
    cdat1 = cdat %>% group_by(group) %>% mutate(adj.Pval=p.adjust(Pval, method="BH"))  
    write.table(cdat1, file=sprintf("%s/HR1/%s/HR.sepCan_%s.txt",dir,i,type), sep="\t", quote=F, col.names=T, row.names=F)
    
    
    cdat2 = read.table(sprintf("%s/HR1/%s/HR.sepCan_%s.txt",dir,i,type), sep="\t", header=T, stringsAsFactors=F)
    cdat2$HR.p = cdat2$HR
    cdat2$HR.p = ifelse(cdat2$HR.p>2, 2, cdat2$HR.p)
    cdat2$HR.range.p = sprintf("%.1f\n%s", cdat2$HR, cdat2$HR.range)
    
    cdat2[!is.na(cdat2$HR) & abs(cdat2$HR)>1000, "HR"] = NA
    cdat2[is.na(cdat2$HR), "HR.range.p"] = "Inf"
    cdat2[cdat2$HR.range.p=="Inf","HR.p"] = 1
    
    cdat2$sig = ifelse(cdat2$adj.Pval <0.05, "adj.P<0.05", NA)
    cdat2$sig = factor(cdat2$sig, levels=c("adj.P<0.05"))
    sig.color = c("#EE2C2C")
    names(sig.color) = c("adj.P<0.05")
    
    p = ggplot(cdat2, aes(x=cancerType, y=group, fill=HR.p)) +
      geom_tile(size=3) +
      geom_tile(aes(color=sig),size=1.2, alpha=0) +
      geom_text(aes(label=HR.range.p), angle=0, color="black", size=3) +
      scale_fill_gradient2( high="#EEAEEE", mid="#FFFFFF", low="#1C86EE", midpoint=1, limits=c(0,2), name="HR", na.value="#FFFFFF") +
      scale_colour_manual(name="significance",values=sig.color) +
      theme_classic() + 
      theme(axis.text.x=element_text(vjust=1, hjust=0, angle=315)) +
      xlab('cancerType') + ylab("group") + ggtitle(sprintf("baseline: %s", levels(dat$immuneType)[1]))
    
    width = 5 + length(unique(cdat2$cancerType))*0.55
    height = 1.75 + length(unique(cdat2$cellType))*0.55
    ggsave(p, file=sprintf("%s/HR1/%s/COX.sepCan_%s.pdf",dir,i,type), width=width, height=height, limitsize=F)
    
    
    ###meta
    ssdir = sprintf("%s/HR1/%s/meta/",dir,i)
    dir.create(ssdir, F, T)
    
    cdat3 = read.table(sprintf("%s/HR1/%s/HR.sepCan_%s.txt",dir,i,type), sep="\t", header=T, stringsAsFactors=F)
    cdat3 = cdat3[cdat3$cancerType!='panCan',]
    cdat3 = cdat3[cdat3$cancerType!='TGCT',]
    cdat3 = cdat3[order(cdat3[,3]),]
    
    meta = meta::metagen( TE=cdat3$coef, seTE=cdat3$coef.se, studlab=cdat3$cancerType,
                          comb.fixed=F, comb.random=T, prediction=F, sm="HR")
    pdf(sprintf("%s/%s.pdf",ssdir,type), width=14, height=10)
    meta::forest(meta, layout="JAMA", test.overall.random=T, digits.pval=4,
                 colgap.forest.left="5cm", zero.pval=T)
    dev.off()
  }
  
}

####plot HR

datas = read.table('HR.sepCan_reg.txt',sep="\t", header=T, stringsAsFactors=F)

datas$HR_low = as.numeric(sapply(strsplit(gsub("[()]", "", datas$HR.range), "-"), "[[", 1))
datas$HR_high = as.numeric(sapply(strsplit(gsub("[()]", "", datas$HR.range), "-"), "[[", 2))
datas$stars <- paste0(
  ifelse(datas$Pval < 0.05, "*", ""),
  ifelse(datas$Pval < 0.01, "*", ""),
  ifelse(datas$Pval < 0.001, "*", "")
)

datas <- datas %>% dplyr::arrange(HR)

rownames(datas)=datas$cancerType
datas=datas[c('SKCM','UVM','KIRC','CRC','LUAD','BLCA','HNSC','LIHC'
              ,'STAD','PAAD','MESO', 'BRCA','ESCA','panCan'),]


datas$sig <- ifelse(datas$Pval < 0.05, "sig","non-sig")
col <- c("sig"="#DB3B3B", "non-sig"="#34326B")

####Plots
rangeb <- range(datas$HR_low, datas$HR_high, na.rm = TRUE)
breaks <- axisTicks(rangeb / 2, log = F, nint = 7)
rangeplot <- rangeb
# make plot twice as wide as needed to create space for annotations
rangeplot[1] <- rangeplot[1] - diff(rangeb)
# increase white space on right for p-vals:
rangeplot[2] <- rangeplot[2] + .15 * diff(rangeb)

y_stars <- rangeb[2]
x_annotate <- seq_len(nrow(datas))

annot_size_mm <- 0.7 * as.numeric(grid::convertX(unit(theme_get()$text$size, "pt"), "mm"))

datas$cancerType <- factor(datas$cancerType)
dev.off()
p <- ggplot(datas, aes(cancerType, HR, color=sig)) +
  geom_point(pch = 15, size = 4) +
  geom_errorbar(aes(ymin = HR_low, ymax = HR_high), width = 0.15) +
  geom_hline(yintercept = 1, linetype = 3) +
  scale_color_manual(values = col) + 
  scale_y_log10(
    name = "Harzard Ratio (95% CI)",
    labels = sprintf("%g", breaks),
    expand = c(0.02, 0.02),
    breaks = breaks
  ) +
  theme_light() +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  ) +
  annotate(
    geom = "text", x = x_annotate, y = y_stars,
    label = datas$stars, size = annot_size_mm,
    fontface = "italic"
  ) 
ggsave("HR.pdf", p, width = 5, height = 3)