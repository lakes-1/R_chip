if (T) {
  rm(list = ls()) 
  options(stringsAsFactors = F)
  load(file = 'data/deg.Rdata')
  head(deg)
}
# 数据预处理
## 不同的阈值，筛选到的差异基因数量就不一样，后面的超几何分布检验结果就大相径庭。
logFC_t=1.5
# 预处理1
if (T) {
  deg$g=ifelse(deg$P.Value>0.05,'stable',
               ifelse( deg$logFC > logFC_t,'UP',
                       ifelse( deg$logFC < -logFC_t,'DOWN','stable') )
  )
  table(deg$g)
  head(deg)
  deg$symbol=rownames(deg)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  df <- bitr(unique(deg$symbol), fromType = "SYMBOL",
             toType = c( "ENTREZID"),
             OrgDb = org.Hs.eg.db)
  head(df)
  DEG=deg
  head(DEG)
  
  DEG=merge(DEG,df,by.y='SYMBOL',by.x='symbol')
  head(DEG)
  save(DEG,file = 'data/anno_DEG.Rdata')
}
# 预处理2
if (T) {
  gene_up= DEG[DEG$g == 'UP','ENTREZID'] 
  gene_down=DEG[DEG$g == 'DOWN','ENTREZID'] 
  gene_diff=c(gene_up,gene_down)
  gene_all=as.character(DEG[ ,'ENTREZID'] )
  data(geneList, package="DOSE") 
  head(geneList)
  boxplot(geneList)
  boxplot(DEG$logFC)
  
  geneList=DEG$logFC
  names(geneList)=DEG$ENTREZID
  geneList=sort(geneList,decreasing = T)
}

# detailed plot,太过于细节，但是部分出图损失信息
if (F) {
  source('kegg_and_go_up_and_down.R')
  run_kegg(gene_up,gene_down,pro='npc_VS_normal')
  # 需要多go数据库的3个条目进行3次富集分析，非常耗时。
  run_go(gene_up,gene_down,pro='npc_VS_normal')
}


# 综合显示图
if (T) {
  go <- enrichGO(gene_up, OrgDb = "org.Hs.eg.db", ont="all") 
  library(ggplot2)
  library(stringr)
  barplot(go, split="ONTOLOGY")+ facet_grid(ONTOLOGY~., scale="free") 
  barplot(go, split="ONTOLOGY",font.size =10)+ 
    facet_grid(ONTOLOGY~., scale="free") + 
    scale_x_discrete(labels=function(x) str_wrap(x, width=50))+
    ggsave('pic/gene_up_GO_all_barplot.png') 
  
  go <- enrichGO(gene_down, OrgDb = "org.Hs.eg.db", ont="all") 
  barplot(go, split="ONTOLOGY",font.size =10)+ 
    facet_grid(ONTOLOGY~., scale="free") + 
    scale_x_discrete(labels=function(x) str_wrap(x, width=50))+
    ggsave('pic/gene_down_GO_all_barplot.png')
}
