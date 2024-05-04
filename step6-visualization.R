
rm(list = ls()) 
options(stringsAsFactors = F)

library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)

load(file = 'anno_DEG.Rdata') 

gene_up= DEG[DEG$g == 'UP','ENTREZID'] 
gene_down=DEG[DEG$g == 'DOWN','ENTREZID'] 
gene_diff=c(gene_up,gene_down)
gene_all=as.character(DEG[ ,'ENTREZID'] ) 
gene_down
gene_up

enrichKK <- enrichKEGG(gene         =  gene_up,
                      organism     = 'hsa',
                      #universe     = gene_all,
                      pvalueCutoff = 0.1,
                      qvalueCutoff =0.1)
head(enrichKK)[,1:6] 
browseKEGG(enrichKK, 'hsa04512')
dotplot(enrichKK)
enrichKK=DOSE::setReadable(enrichKK, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
enrichKK 

#(3)visualization

par(mfrow=c(2,1))
barplot(enrichKK,showCategory=20)

dotplot(enrichKK)

geneList = deg$logFC
names(geneList)=deg$ENTREZID
geneList = sort(geneList,decreasing = T)

#Gene-Concept Network
cnetplot(enrichKK, categorySize="pvalue", foldChange=geneList,colorEdge = TRUE)
cnetplot(enrichKK, foldChange=geneList, circular = TRUE, colorEdge = TRUE)
#Enrichment Map
emapplot(enrichKK)

# goplot(enrichKK)
#(5)Heatmap-like functional classification
heatplot(enrichKK,foldChange = geneList)


