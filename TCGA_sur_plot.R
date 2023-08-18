library(AUCell)
library(Seurat)
library(SingleCellExperiment)
suppressPackageStartupMessages({
  library("plyr")
  library("dplyr")
  library("ggpubr")
  library("ggsci")
  library("reshape2")
  library("survival")
  library("survminer")
  library("data.table")
})

sce = readRDS('./TCGA_Toil.rds')
sce = sce[,sce$sampleType=="Primary.Solid.Tumor"]
se = as.Seurat(sce,data='exprs')
geneSets = list(TaNK=c('KLRF1', 'DNAJB1', 'BAG3', 'SERPINE1' ,'NR4A1'))
new = 'True'

for (eachc in unique(se$cancer.type.abbreviation)) {
  se.c = subset(se,cancer.type.abbreviation==eachc)
  
  cells_rankings <- AUCell_buildRankings(se.c@assays$originalexp@data)
  
  cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
  aucs <- as.numeric(getAUC(cells_AUC)['TaNK', ])
  
  
  bdf = data.frame(sample=se.c$sample,cancer.type.abbreviation=eachc,
                   aucs=aucs)
  if (new=='True') {
    adf = bdf
    new='added'
  }else{
    adf = rbind(adf,bdf)
  }
}


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

dir="0.out"

###KM

dir.create(sprintf("%s/KM/",dir), F, T)
  
for (cancer in c("panCan", unique(dat$cancerType))){
  pdat = dat
  if(cancer!="panCan"){
    pdat = dat[dat$cancerType==cancer,]
  }
  res.cut <- surv_cutpoint(pdat, time = "OS.time.use", 
                             event = "OS.use", 
                             variables = i)
  res.cat <- surv_categorize(res.cut)
  pdat$group <- res.cat[,i] 
  pdat$group = factor(pdat$group)
    
  fit = survfit(Surv(OS.time.use,OS.use) ~ group, data=pdat)
    
  p = ggsurvplot(fit, pdat, size=0.3, vlegend.labs=unique(pdat$group),
                   surv.median.line="none", pval=T, conf.int=F,
                   #risk.table=T, risk.table.y.text.col=T,
                   palette=c("#990066","#CCCCCC"),title=cancer) + 
                   xlab("Years")

  pdf(file=sprintf("%s/KM/KM.%s.pdf", dir, cancer), width=3, height=3,onefile = FALSE)
  print(p)
  dev.off()
}
