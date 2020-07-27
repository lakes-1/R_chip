rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)
# 下面代码只需要修改file的名字
# 下载的文件会自动保存到./data/下
# 得到的gset是一个ExpressionSet对象
# 在注释的过程中需注意数据对应的平台

# 只需要修改file的名字，即可下载得到相应的geo数据
# 获取患者信息，需要自行修改
file='GSE64634'
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64634

if (T) {
  # 数据下载
  if (T) {
    library(GEOquery)
    # 这个包需要注意两个配置，一般来说自动化的配置是足够的。
    #Setting options('download.file.method.GEOquery'='auto')
    #Setting options('GEOquery.inmemory.gpl'=FALSE)
    fdata=paste0(file,"_eSet.Rdata")
    fpath=paste0("data/",fdata)
    if(!file.exists(fpath)){
      gset <- getGEO(file, destdir="data/",
                     AnnotGPL = F,     ## 注释文件
                     getGPL = F)       ## 平台文件
      gset=gset[[1]]
      save(gset,file=fpath)   ## 保存到本地
    }
    load(fpath)
  }
  gset
  
  # 获取患者信息，这里需要自行修改
  if (T) {
    pd=pData(gset)# 根据disease state:ch1一列得知分组信息
    group_list=c(rep('normal',4),rep('npc',12))# nasopharyngeal carcinoma NPC 鼻咽癌
    table(group_list)
  }
  
  # 对数据进行normalization，用的是limma的片内标准化
  if (T) {
    dat=exprs(gset)
    dim(dat)
    
    dat[1:4,1:4]
    boxplot(dat,las=2)
    dat=dat[apply(dat,1,sd)>0,]# 去除都是0的探针
    dat[dat<0]=1
    boxplot(dat,las=2)
    
    #dat=log2(dat+1)
    #boxplot(dat,las=2)
    library(limma)
    dat=normalizeBetweenArrays(dat)
    boxplot(dat,las=2)
  }
  
  
  # 探针注释2种方法，推荐方法2，但我觉得可以用提取到的芯片平台来做批量化的处理，下载其注释信息
  # 方法1:比较麻烦，而且不方便，一般不用这种方法
  if(F){
    library(GEOquery)
    #Download GPL file, put it in the current directory, and load it:
    gpl <- getGEO('GPL570', destdir="data/")
    GPL=Table(gpl)
    colnames(GPL)
    head(GPL[,c(1,11)]) ## you need to check this , which column do you need
    probe2gene=GPL[,c(1,11)]
    head(probe2gene)
    save(probe2gene,file='probe2gene.Rdata')
  }
  
  # 方法2:用hgu133plus2.db这个R包比较方便,这是已知知道了GPL平台且知道了平台对应的包所选用的注释信息
  if (T) {
    library(hgu133plus2.db)
    ids=toTable(hgu133plus2SYMBOL)
    head(ids)
    ids=ids[ids$symbol != '',]
    ids=ids[ids$probe_id %in%  rownames(dat),]# 过滤没法注释的探针
    
    dat[1:4,1:4]   
    dat=dat[ids$probe_id,]# 调整顺序，让dat的顺序和ids中的一致
    
    ids$median=apply(dat,1,median)
    ids=ids[order(ids$symbol,ids$median,decreasing = T),]# 按照基因名、中位数大小排序
    ids=ids[!duplicated(ids$symbol),]# 只保留相同symbol中中位数最大的探针
    dat=dat[ids$probe_id,]# 调整顺序，让dat的顺序和ids中的一致
    rownames(dat)=ids$symbol# id转换
    dat[1:4,1:4]
  }
  
  # 方法3：用jimmy的idmap系列包，1版本是针对bioconductor,2版本是针对大多数GPL平台的，
  #3版本是针对一些以探针序列信息的包
  #其中一般来说只用于一个gse号只有一个平台，但有多个平台也可以这么做，循环即可
  #调用ls('package:idmap2')，idmap2直接观察平台即可，idmap1则是通过p2s_df的提取
  if (T) {
    platform=gset[[1]]@annotation
    ids=p2s_df[p2s_df[,3]%in%platform,]
    ids=ids[ids$symbol != '',]
    ids=ids[ids$probe_id %in%  rownames(dat),]# 过滤没法注释的探针
    dat[1:4,1:4]   
    dat=dat[ids$probe_id,]# 调整顺序，让dat的顺序和ids中的一致
    ids$median=apply(dat,1,median)
    ids=ids[order(ids$symbol,ids$median,decreasing = T),]# 按照基因名、中位数大小排序
    ids=ids[!duplicated(ids$symbol),]# 只保留相同symbol中中位数最大的探针
    dat=dat[ids$probe_id,]# 调整顺序，让dat的顺序和ids中的一致
    rownames(dat)=ids$symbol# id转换
    dat[1:4,1:4]
  }
save(dat,group_list,file = 'data/step1-output.Rdata')