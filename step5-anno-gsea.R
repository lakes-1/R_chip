# 载入数据和R包
# DEG为limma得到的差异分析结果
if (T) {
  rm(list = ls()) 
  options(stringsAsFactors = F)
  load(file = 'data/anno_DEG.Rdata')
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
}


### 对 MigDB中的全部基因集 做GSEA分析。
### 按照FC的值对差异基因进行排序
# http://www.bio-info-trainee.com/2105.html
# http://www.bio-info-trainee.com/2102.html 
# 自行修改存放gmt文件路径d
# GSEA每个gene set的具体结果保存在gsea_results这个list中
# 而最终结果保存在gsea_results_df数据框中
d='D:/搜狗高速下载/msigdb.v7.0.symbols/'
if(T){
  geneList=DEG$logFC
  names(geneList)=DEG$symbol
  geneList=sort(geneList,decreasing = T)
  #选择gmt文件（MigDB中的全部基因集）
  
  #gmts=list.files(d,pattern = 'all')
  gmts=list.files(d)
  gmts
  
  #GSEA分析
  library(GSEABase) # BiocManager::install('GSEABase')
  ## 下面使用lapply循环读取每个gmt文件，并且进行GSEA分析
  ## 如果存在之前分析后保存的结果文件，就不需要重复进行GSEA分析。
  f='data/gsea_results.Rdata'
if(!file.exists(f)){
    gsea_results <- lapply(gmts, function(gmtfile){
      # gmtfile=gmts[2]

      filepath=paste0(d,gmtfile)
      geneset <- read.gmt(filepath) 
      print(paste0('Now process the ',gmtfile))
      egmt <- GSEA(geneList, TERM2GENE=geneset, verbose=FALSE)
      head(egmt)
      # gseaplot(egmt, geneSetID = rownames(egmt[1,]))
      return(egmt)
    })
    # 上面的代码耗时，所以保存结果到本地文件
    save(gsea_results,file = f)
  }
  load(file = f)
  #提取gsea结果，熟悉这个对象
  gsea_results_list<- lapply(gsea_results, function(x){
    cat(paste(dim(x@result)),'\n')
    x@result
  })
}
# 随便看几个结果图,canshukeyijiatitle
dev.off()
gsea_results_df <- do.call(rbind, gsea_results_list)
gseaplot(gsea_results[[1]],geneSetID = "NIKOLSKY_BREAST_CANCER_7P15_AMPLICON") 
gseaplot(gsea_results[[1]],'RICKMAN_HEAD_AND_NECK_CANCER_D',) 
#