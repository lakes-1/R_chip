### 主要是关于KEGG方面的扩展图
### 主要是关于KEGG方面的扩展图
# 载入数据
if (T) {
  rm(list = ls()) 
  options(stringsAsFactors = F)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  load(file = 'data/anno_DEG.Rdata')  
  
  head(DEG)
}

# pre-process data
if (T) {
  gene_up= DEG[DEG$g == 'UP','ENTREZID'] 
  gene_down=DEG[DEG$g == 'DOWN','ENTREZID'] 
  gene_diff=c(gene_up,gene_down)
  gene_all=as.character(DEG[ ,'ENTREZID'] )
  # 制作差异基因list L
  boxplot(DEG$logFC)
  
  geneList=DEG$logFC
  names(geneList)=DEG$ENTREZID
  geneList=sort(geneList,decreasing = T)
}


# KEGG富集分析得到结果,这边只做了上调
if (T) {
  if (!file.exists("data/enrichkk.rdata")) {
    gene_down
    gene_up
    enrichKK <- enrichKEGG(gene         =  gene_up,
                           organism     = 'hsa',
                           #universe     = gene_all,
                           pvalueCutoff = 0.1,
                           qvalueCutoff =0.1)
    save(enrichKK,file = "data/enrichkk.rdata")
  }
  load(file = "data/enrichkk.rdata")
  # 查看KEGG结果
  head(enrichKK)[,1:6] 
  # 打开网页看相关KEGG通路图
  browseKEGG(enrichKK, 'hsa05150')
  
  # 将数据中的entrz-id变成symbol
  # 更为易读
  enrichKK=DOSE::setReadable(enrichKK, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
  enrichKK 
}


## 可视化
#条带图
if (T) {
  # par(mfrow=c(2,1))
  barplot(enrichKK,showCategory=20)
  ggsave("pic/barplot.png")
}

#气泡图
if (T) {
  dotplot(enrichKK)
  ggsave("pic/dotplot.png")
}

#下面的图需要映射颜色，设置和示例数据一样的geneList

# 展示top5通路的共同基因，要放大看。
#Gene-Concept Network
if (T) {
  cnetplot(enrichKK, foldChange=geneList,colorEdge = TRUE, circular = F)
  ggsave("pic/cnetplot.png")
  cnetplot(enrichKK, foldChange=geneList, colorEdge = TRUE, circular = T)
  ggsave("pic/cnetplot_circular.png")
}


#Enrichment Map
if (T) {
  emapplot(enrichKK)
  ggsave("pic/Enrichment_Map.png")
}

#(4)展示通路关系,仅仅是针对于GO数据库结果。
# goplot(enrichKK)
#(5)Heatmap-like functional classification
if (T) {
  heatplot(enrichKK,foldChange = geneList)
  ggsave("pic/Enrichment_Heatmap.png")
}
