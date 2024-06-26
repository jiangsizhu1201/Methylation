
rm(list = ls())   
#options()$repos 
#options()$BioC_mirror
#options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
#options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
#options()$repos 
#options()$BioC_mirror


# https://bioconductor.org/packages/release/bioc/html/GEOquery.html
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("KEGG.db",ask = F,update = F)

BiocManager::install("minfi",ask = F,update = F)
BiocManager::install("ChAMP",ask = F,update = F)
BiocManager::install("methylationArrayAnalysis",ask = F,update = F) 
BiocManager::install("wateRmelon",ask = F,update = F) 


BiocManager::install(c("GSEABase","GSVA","clusterProfiler" ),ask = F,update = F)
BiocManager::install(c("GEOquery","limma","impute" ),ask = F,update = F)
BiocManager::install(c("org.Hs.eg.db","hgu133plus2.db" ),ask = F,update = F)


# source("https://bioconductor.org/biocLite.R") 
# library('BiocInstaller') 
# options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") 
# BiocInstaller::biocLite("GEOquery")
# BiocInstaller::biocLite(c("limma"))
# BiocInstaller::biocLite(c("impute"))

options()$repos
install.packages('WGCNA')
install.packages(c("FactoMineR", "factoextra"))
install.packages(c("ggplot2", "pheatmap","ggpubr"))
library("FactoMineR")
library("factoextra")

library(GSEABase)
library(GSVA)
library(clusterProfiler)
library(ggplot2)
library(ggpubr)
library(hgu133plus2.db)
library(limma)
library(org.Hs.eg.db)
library(pheatmap)

