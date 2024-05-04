

rm(list = ls())   
options(stringsAsFactors = F)
library("ChAMP")
library("minfi")
require(GEOquery)
require(Biobase)
load(file = 'step1-output.Rdata')

myLoad  

if(F){
  myNorm <- champ.norm(beta=myLoad$beta,arraytype="450K",cores=5)
  dim(myNorm) 
  pD=myLoad$pd
  save(myNorm,pD,file = 'step2-champ_myNorm.Rdata')
}
load(file = 'step2-champ_myNorm.Rdata')
# 
beta.m=myNorm
group_list=myLoad$pd$Group
dim(beta.m) 
# 
if(T){
  
  
  dat=t(beta.m)
  dat[1:4,1:4] 
  library("FactoMineR")
  library("factoextra")  

  dat.pca <- PCA(dat , graph = FALSE) 
  fviz_pca_ind(dat.pca,
               geom.ind = "point", # show points only (nbut not "text")
               col.ind = group_list, # color by groups
               # palette = c("#00AFBB", "#E7B800"),
               addEllipses = TRUE, # Concentration ellipses
               legend.title = "Groups"
  )
  ggsave('all_samples_PCA.png')
  
  dat=beta.m
  dat[1:4,1:4] 
  cg=names(tail(sort(apply(dat,1,sd)),1000))
  library(pheatmap)
  pheatmap(dat[cg,],show_colnames =F,show_rownames = F) 
  n=t(scale(t(dat[cg,]))) 
  n[n>2]=2 
  n[n< -2]= -2
  n[1:4,1:4]
  pheatmap(n,show_colnames =F,show_rownames = F)
  ac=data.frame(group=group_list)
  rownames(ac)=colnames(n)  
  pheatmap(n,show_colnames =F,show_rownames = F,
           annotation_col=ac,filename = 'heatmap_top1000_sd.png')
  dev.off()
  
  exprSet=beta.m
  pheatmap::pheatmap(cor(exprSet)) 

  colD=data.frame(group_list=group_list)
  rownames(colD)=colnames(exprSet)
  pheatmap::pheatmap(cor(exprSet),
                     annotation_col = colD,
                     show_rownames = F,
                     filename = 'cor_all.png')
  dev.off() 
  exprSet=exprSet[names(sort(apply(exprSet, 1,mad),decreasing = T)[1:500]),]
  dim(exprSet)
  # M=cor(log2(exprSet+1)) 
  M=cor(exprSet)
  pheatmap::pheatmap(M,annotation_col = colD)
  pheatmap::pheatmap(M,
                     show_rownames = F,
                     annotation_col = colD,
                     filename = 'cor_top500.png')
  dev.off() 
  
}


myLoad 
beta.m=myLoad$beta

if(F){
  # wateRmelon 
  library("wateRmelon")
  beta.m=beta.m[rowMeans(beta.m)>0.005,]
  pdf(file="rawBox.pdf")
  boxplot(beta.m,col = "blue",xaxt = "n",outline = F)
  dev.off()
  beta.m = betaqn(beta.m)
  pdf(file="normalBox.pdf")
  boxplot(beta.m,col = "red",xaxt = "n",outline = F)
  dev.off()
  
  # 
  group_list=myLoad$pd$Group
  pdf(file="densityBeanPlot.pdf")
  par(oma=c(2,10,2,2))
  densityBeanPlot(beta.m, sampGroups = group_list)
  dev.off()
  pdf(file="mdsPlot.pdf")
  mdsPlot(beta.m, numPositions = 1000, sampGroups = group_list)
  dev.off()
  

  grset=makeGenomicRatioSetFromMatrix(beta.m,what="Beta")
  M = getM(grset)
  # 因为甲基化芯片是450K或者850K
  dmp <- dmpFinder(M, pheno=group_list, type="categorical")
  dmpDiff=dmp[(dmp$qval<0.05) & (is.na(dmp$qval)==F),]
  dim(dmpDiff)
}



