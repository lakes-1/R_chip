### 对 MigDB中的全部基因集 做GSVA分析。
## 还有ssGSEA, PGSEA
# 载入数据
if(T){
  rm(list = ls()) 
  options(stringsAsFactors = F)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  load(file = 'data/step1-output.Rdata')
  # 每次都要检测数据
  dat[1:4,1:4]  
}

# GSVA分析
# 存放gene set的文件路径需要具体修改
d='D:/搜狗高速下载/msigdb.v7.0.symbols/'
if (T) {
  X=dat
  table(group_list)
  ## Molecular Signatures Database (MSigDb) 
  #gmts=list.files(d,pattern = 'all')
  gmts=list.files(d)
  gmts
  library(GSVA) # BiocManager::install('GSVA')
  if(!file.exists('data/gsva_msigdb.Rdata')){
    es_max <- lapply(gmts, function(gmtfile){ 
      # gmtfile=gmts[8];gmtfile
      geneset <- read.gmt(file.path(d,gmtfile))  
      es.max <- gsva(X, geneset, 
                     mx.diff=FALSE, verbose=FALSE, 
                     parallel.sz=1)
      return(es.max)
    })
    adjPvalueCutoff <- 0.001
    logFCcutoff <- log2(2)
    es_deg <- lapply(es_max, function(es.max){
      # es.max=es_max[[1]]
      table(group_list)
      dim(es.max)
      design <- model.matrix(~0+factor(group_list))
      colnames(design)=levels(factor(group_list))
      rownames(design)=colnames(es.max)
      design
      library(limma)
      # contrast.matrix<-makeContrasts(paste0(unique(group_list),collapse = "-"),levels = design)
      contrast.matrix<-makeContrasts("npc-normal",
                                     levels = design)
      
      contrast.matrix ##这个矩阵声明，我们要把progres.组跟stable进行差异分析比较
      
      deg = function(es.max,design,contrast.matrix){
        ##step1
        fit <- lmFit(es.max,design)
        ##step2
        fit2 <- contrasts.fit(fit, contrast.matrix) 
        ##这一步很重要，大家可以自行看看效果
        
        fit2 <- eBayes(fit2)  ## default no trend !!!
        ##eBayes() with trend=TRUE
        ##step3
        res <- decideTests(fit2, p.value=adjPvalueCutoff)
        summary(res)
        tempOutput = topTable(fit2, coef=1, n=Inf)
        nrDEG = na.omit(tempOutput) 
        #write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
        head(nrDEG)
        return(nrDEG)
      }
      
      re = deg(es.max,design,contrast.matrix)
      nrDEG=re
      head(nrDEG) 
      return(nrDEG)
    })
    gmts
    save(es_max,es_deg,file='data/gsva_msigdb.Rdata')
  }
}

# 画图展示，结果存放在pic/下
if (T) {
  load(file='data/gsva_msigdb.Rdata')
  library(pheatmap)
  lapply(1:length(es_deg), function(i){
    # i=8
    print(i)
    dat=es_max[[i]]
    df=es_deg[[i]]
    df=df[df$P.Value<0.01 & abs(df$logFC) > 0.3,]
    print(dim(df))
    if(nrow(df)>5){
      n=rownames(df)
      dat=dat[match(n,rownames(dat)),]
      ac=data.frame(g=group_list)
      rownames(ac)=colnames(dat)
      rownames(dat)=substring(rownames(dat),1,50)
      pheatmap::pheatmap(dat, 
                         fontsize_row = 8,height = 11,
                         annotation_col = ac,show_colnames = F,
                         filename = paste0('[pic/gsva_',strsplit(gmts[i],'[.]')[[1]][1],'.pdf'))
      
    }
  })
  
  adjPvalueCutoff <- 0.001
  logFCcutoff <- log2(2)
  df=do.call(rbind ,es_deg)
  es_matrix=do.call(rbind ,es_max)
  df=df[df$P.Value<0.01 & abs(df$logFC) > 0.5,]
  write.csv(df,file = 'data/GSVA_DEG.csv')
}
