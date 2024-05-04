## KEGG pathway analysis
### 
run_kegg <- function(gene_up,gene_down,geneList=F,pro='test'){
  gene_up=unique(gene_up)
  gene_down=unique(gene_down)
  gene_diff=unique(c(gene_up,gene_down))
  ###   over-representation test

  kk.up <- enrichKEGG(gene         = gene_up,
                      organism     = 'hsa',
                      #universe     = gene_all,
                      pvalueCutoff = 0.9,
                      qvalueCutoff =0.9)
  head(kk.up)[,1:6]
  kk=kk.up
  dotplot(kk)
  kk=DOSE::setReadable(kk, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
  write.csv(kk@result,paste0(pro,'_kk.up.csv'))
  

  kk.down <- enrichKEGG(gene         =  gene_down,
                        organism     = 'hsa',
                        #universe     = gene_all,
                        pvalueCutoff = 0.9,
                        qvalueCutoff =0.9)
  head(kk.down)[,1:6]
  kk=kk.down
  dotplot(kk)
  kk=DOSE::setReadable(kk, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
  write.csv(kk@result,paste0(pro,'_kk.down.csv'))
  

  kk.diff <- enrichKEGG(gene         = gene_diff,
                        organism     = 'hsa',
                        pvalueCutoff = 0.05)
  head(kk.diff)[,1:6]
  kk=kk.diff
  dotplot(kk)
  kk=DOSE::setReadable(kk, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
  write.csv(kk@result,paste0(pro,'_kk.diff.csv'))
  
  
  kegg_diff_dt <- as.data.frame(kk.diff)
  kegg_down_dt <- as.data.frame(kk.down)
  kegg_up_dt <- as.data.frame(kk.up)
  down_kegg<-kegg_down_dt[kegg_down_dt$pvalue<0.01,];down_kegg$group=-1
  up_kegg<-kegg_up_dt[kegg_up_dt$pvalue<0.01,];up_kegg$group=1

  g_kegg=kegg_plot(up_kegg,down_kegg)
  print(g_kegg)
  
  ggsave(g_kegg,filename = paste0(pro,'_kegg_up_down.png') )
  
if(geneList){

  kk_gse <- gseKEGG(geneList     = geneList,
                    organism     = 'hsa',
                    nPerm        = 1000,
                    minGSSize    = 20,
                    pvalueCutoff = 0.9,
                    verbose      = FALSE)
  head(kk_gse)[,1:6]
  gseaplot(kk_gse, geneSetID = rownames(kk_gse[1,]))
  gseaplot(kk_gse, 'hsa04110',title = 'Cell cycle') 
  kk=DOSE::setReadable(kk_gse, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
  tmp=kk@result
  write.csv(kk@result,paste0(pro,'_kegg.gsea.csv'))
  
  
  down_kegg<-kk_gse[kk_gse$pvalue<0.05 & kk_gse$enrichmentScore < 0,];down_kegg$group=-1
  up_kegg<-kk_gse[kk_gse$pvalue<0.05 & kk_gse$enrichmentScore > 0,];up_kegg$group=1
  
  g_kegg=kegg_plot(up_kegg,down_kegg)
  print(g_kegg)
  ggsave(g_kegg,filename = paste0(pro,'_kegg_gsea.png'))
  # 
}
  
}

### GO database analysis 

run_go <- function(gene_up,gene_down,pro='test'){
  gene_up=unique(gene_up)
  gene_down=unique(gene_down)
  gene_diff=unique(c(gene_up,gene_down))
  g_list=list(gene_up=gene_up,
              gene_down=gene_down,
              gene_diff=gene_diff)

  if(T){
    go_enrich_results <- lapply( g_list , function(gene) {
      lapply( c('BP','MF','CC') , function(ont) {
        cat(paste('Now process ',ont ))
        ego <- enrichGO(gene          = gene,
                        #universe      = gene_all,
                        OrgDb         = org.Hs.eg.db,
                        ont           = ont ,
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.99,
                        qvalueCutoff  = 0.99,
                        readable      = TRUE)
        
        print( head(ego) )
        return(ego)
      })
    })
    save(go_enrich_results,file =paste0(pro, '_go_enrich_results.Rdata'))
    
  }

  load(file=paste0(pro, '_go_enrich_results.Rdata'))
  
  n1= c('gene_up','gene_down','gene_diff')
  n2= c('BP','MF','CC') 
  for (i in 1:3){
    for (j in 1:3){
      fn=paste0(pro, '_dotplot_',n1[i],'_',n2[j],'.png')
      cat(paste0(fn,'\n'))
      png(fn,res=150,width = 1080)
      print( dotplot(go_enrich_results[[i]][[j]] ))
      dev.off()
    }
  }
  
  
}

kegg_plot <- function(up_kegg,down_kegg){
  dat=rbind(up_kegg,down_kegg)
  colnames(dat)
  dat$pvalue = -log10(dat$pvalue)
  dat$pvalue=dat$pvalue*dat$group 
  
  dat=dat[order(dat$pvalue,decreasing = F),]
  
  g_kegg<- ggplot(dat, aes(x=reorder(Description,order(pvalue, decreasing = F)), y=pvalue, fill=group)) + 
    geom_bar(stat="identity") + 
    scale_fill_gradient(low="blue",high="red",guide = FALSE) + 
    scale_x_discrete(name ="Pathway names") +
    scale_y_continuous(name ="log10P-value") +
    coord_flip() + theme_bw()+theme(plot.title = element_text(hjust = 0.5))+
    ggtitle("Pathway Enrichment") 
}


