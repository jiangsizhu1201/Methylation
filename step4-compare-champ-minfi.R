
rm(list = ls())   
options(stringsAsFactors = F)
library("ChAMP")
library("minfi")
require(GEOquery)
require(Biobase)
load(file = 'step1-output.Rdata')
beta.m=myLoad$beta

grset=makeGenomicRatioSetFromMatrix(beta.m,what="Beta")
M = getM(grset)
group_list=myLoad$pd$Group

dmp <- dmpFinder(M, pheno=group_list, type="categorical")
dmpDiff=dmp[(dmp$qval<0.05) & (is.na(dmp$qval)==F),]
dim(dmpDiff)

load(file = 'step3-output-myDMP.Rdata')
champDiff=myDMP[[1]]

dim(dmpDiff)
dim(champDiff)
length(intersect(rownames(dmpDiff),rownames(champDiff)))

source('functions_DEM.R')
visual_champ_DEM(myLoad,myDMP,group_list,pro='test')



