

rm(list = ls())   
options(stringsAsFactors = F)
library("ChAMP")
library("minfi")
require(GEOquery)
require(Biobase)
load(file = 'step1-output.Rdata')
myLoad    # 
load(file = 'step2-champ_myNorm.Rdata')
group_list=myLoad$pd$Group
table(group_list)
myDMP <- champ.DMP(beta = myNorm,pheno=group_list)
head(myDMP[[1]])
save(myDMP,file = 'step3-output-myDMP.Rdata')

# DMP.GUI(DMP=myDMP[[1]],beta=myNorm,group_list)


if(F){
  myDMR <- champ.DMR(beta = myNorm,pheno=group_list,method="Bumphunter")
  DMR.GUI(DMR=myDMR)
  
  myBlock <- champ.Block(beta = myNorm,pheno=group_list,arraytype="450K")
  head(myBlock$Block)
  Block.GUI(Block=myBlock,beta = myNorm,pheno=group_list,
            runDMP=TRUE,compare.group=NULL,arraytype="450K")
  
  myGSEA <- champ.GSEA(beta=myNorm,DMP=myDMP[[1]],
                       DMR=myDMR, arraytype="450K",adjPval=0.05, method="fisher")
  
  head(myGSEA$DMP)
  head(myGSEA$DMR)
  myEpiMod <- champ.EpiMod(beta=myNorm,pheno=group_list)
  
}




