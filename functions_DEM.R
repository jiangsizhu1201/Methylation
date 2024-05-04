visual_champ_DEM <- function(myLoad,myDMP,group_list,pro='test'){
  beta.m=myLoad$beta
  champDiff=myDMP[[1]] 
  head(champDiff)  
  colnames(champDiff)
  ## for volcano 
  if(T){
    nrDEG=champDiff
    head(nrDEG)
    attach(nrDEG)
    plot(logFC,-log10(P.Value))
    library(ggpubr)
    df=nrDEG
    df$v= -log10(P.Value) #df新增加一列'v',值为-log10(P.Value)
    ggscatter(df, x = "logFC", y = "v",size=0.5)
    
    df$g=ifelse(df$P.Value>0.05,'stable', #if P.Value>0.01，stable gene
                ifelse( df$logFC >0.5,'up', 
                        ifelse( df$logFC < -0.5,'down','stable') )
    )
    table(df$g)
    df$name=rownames(df)
    head(df)
    ggscatter(df, x = "logFC", y = "v",size=0.5,color = 'g')
    ggscatter(df, x = "logFC", y = "v", color = "g",size = 0.5,
              label = "name", repel = T,
              #label.select = rownames(df)[df$g != 'stable'] ,
              label.select =  head( rownames(df)[df$g != 'stable']),
              palette = c("#00AFBB", "#E7B800", "#FC4E07") )
    ggsave(paste0(pro,'_volcano.png'))
    
    ggscatter(df, x = "AveExpr", y = "logFC",size = 0.2)
    df$p_c = ifelse(df$P.Value<0.001,'p<0.001',
                    ifelse(df$P.Value<0.01,'0.001<p<0.01','p>0.01'))
    table(df$p_c )
    ggscatter(df,x = "AveExpr", y = "logFC", color = "p_c",size=0.2, 
              palette = c("green", "red", "black") )
    ggsave(paste0(pro,'_MA.png'))
    
    
  }
  
  ## for heatmap 
  if(T){  
    dat=beta.m
    dat[1:4,1:4]
    table(group_list)
    deg=champDiff
    x=deg$logFC 
    names(x)=rownames(deg) 
    cg=c(names(head(sort(x),100)),
         names(tail(sort(x),100)))
    library(pheatmap)
    pheatmap(dat[cg,],show_colnames =F,show_rownames = F) 
    n=t(scale(t(dat[cg,])))
    
    n[n>2]=2
    n[n< -2]= -2
    n[1:4,1:4]
    pheatmap(n,show_colnames =F,show_rownames = F)
    ac=data.frame(group=group_list)
    rownames(ac)=colnames(n) 
    pheatmap(n,show_colnames =F,
             show_rownames = F,
             cluster_cols = T, 
             annotation_col=ac,
             filename = paste0(pro,'_heatmap_top200_DEG_scale.png'))  
    pheatmap(dat[cg,],show_colnames =F,
             show_rownames = F,
             cluster_cols = T, 
             annotation_col=ac,
             filename = paste0(pro,'_heatmap_top200_DEG_raw.png'))  
    
  }
  
 
}
