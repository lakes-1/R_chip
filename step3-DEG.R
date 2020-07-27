# 载入数据,检查数据
if (T) {
  rm(list = ls()) 
  options(stringsAsFactors = F)
  library(ggpubr)
  load(file = 'data/step1-output.Rdata')
  # 每次都要检测数据
  dat[1:4,1:4] 
  table(group_list)
  boxplot(dat[1,]~group_list) #按照group_list分组画箱线图
  
  # boxplot的美化版
  bplot=function(g){
    df=data.frame(gene=g,stage=group_list)
    p <- ggboxplot(df, x = "stage", y = "gene",
                   color = "stage", palette = "jco",
                   add = "jitter")
    #  Add p-value
    p + stat_compare_means()
  }
}
# 利用定义好的函数检查数据
bplot(dat[1,])
bplot(dat[2,])
bplot(dat[3,])
bplot(dat[4,])


# limma
library(limma)
# 方法1：不制作比较矩阵，简单
# 但是做不到随心所欲的指定任意两组进行比较
if (T) {
  design=model.matrix(~factor( group_list ))
  fit=lmFit(dat,design)
  fit2=eBayes(fit)
  ## 上面是limma包用法的一种方式 
  options(digits = 4) #设置全局的数字有效位数为4
  topTable(fit2,coef=2,adjust='BH') 
}

# 方法2：制作比较矩阵
# 可以随心所欲的指定任意两组进行比较
if (T) {
  design <- model.matrix(~0+factor(group_list))
  colnames(design)=levels(factor(group_list))
  head(design)
  exprSet=dat
  rownames(design)=colnames(exprSet)
  design
  
  # 比较矩阵
  # 这个矩阵声明，我们要把 npc 组跟 Normal 进行差异分析比较
  contrast.matrix<-makeContrasts("npc-normal",
                                 levels = design)
  contrast.matrix
  
  deg = function(exprSet,design,contrast.matrix){
    ##step1
    fit <- lmFit(exprSet,design)
    ##step2
    fit2 <- contrasts.fit(fit, contrast.matrix)
    
    fit2 <- eBayes(fit2)  ## default no trend !!!
    ##eBayes() with trend=TRUE
    
    ##step3
    # 有了比较矩阵后，coef=1，而number=Inf是把所有结果都打印出来
    tempOutput = topTable(fit2, coef=1, number =Inf)
    nrDEG = na.omit(tempOutput) 
    #write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
    head(nrDEG)
    return(nrDEG)
  }
  
  deg = deg(exprSet,design,contrast.matrix)
  
  head(deg)
}


save(deg,file = 'data/deg.Rdata')

load(file = 'data/deg.Rdata')
head(deg)
bplot(dat[rownames(deg)[1],])
## for volcano and MA plot
# 结果存放在pic/volcano.png和pic/MA.png
if(T){
  nrDEG=deg
  head(nrDEG)
  attach(nrDEG)
  # 原始版火山图
  plot(logFC,-log10(P.Value))
  library(ggpubr)
  df=nrDEG
  df$y= -log10(P.Value)
  ggscatter(df, x = "logFC", y = "y",size=0.5)
  # 定义logFC=2为阈值
  df$state=ifelse(df$P.Value>0.01,'stable',
                  ifelse( df$logFC >2,'up',
                          ifelse( df$logFC < -2,'down','stable') )
  )
  table(df$state)
  df$name=rownames(df)
  head(df)
  ggscatter(df, x = "logFC", y = "y",size=0.5,color = 'state')
  ggscatter(df, x = "logFC", y = "y", color = "state",size = 0.5,
            label = "name", repel = T,
            #label.select = rownames(df)[df$state != 'stable'] ,
            label.select = c('TTC9', 'AQP3', 'CXCL11','PTGS2'), #挑选一些基因在图中显示出来
            palette = c("#00AFBB", "#E7B800", "#FC4E07"))
  ggsave('pic/volcano.png')
  
  # MA图
  ggscatter(df, x = "AveExpr", y = "logFC",size = 0.2)
  df$p_c = ifelse(df$P.Value<0.001,'p<0.001',
                  ifelse(df$P.Value<0.01,'0.001<p<0.01','p>0.01'))
  table(df$p_c )
  ggscatter(df,x = "AveExpr", y = "logFC", color = "p_c",size=0.2, 
            palette = c("green", "red", "black") )
  ggsave('pic/MA.png')
}

## for heatmap 
if(T){ 
  load(file = 'data/step1-output.Rdata')
  # 每次都要检测数据
  dat[1:4,1:4]
  table(group_list)
  x=deg$logFC
  names(x)=rownames(deg)
  
  # cg中存放着变化上升和下降的前100个基因名
  cg=c(names(head(sort(x),100)),
       names(tail(sort(x),100)))
  library(pheatmap)
  n=t(scale(t(dat[cg,])))
  
  n[n>2]=2
  n[n< -2]= -2
  n[1:4,1:4]
  pheatmap(n,show_colnames =F,show_rownames = F)
  ac=data.frame(group=group_list)
  rownames(ac)=colnames(n) #将ac的行名也就分组信息（是‘no TNBC’还是‘TNBC’）给到n的列名，即热图中位于上方的分组信息
  pheatmap(n,show_colnames =F,
           show_rownames = F,
           cluster_cols = F, 
           annotation_col=ac,filename = 'pic/heatmap_top200_DEG.png') #列名注释信息为ac即分组信息
  
  
}

write.csv(deg,file = 'data/deg.csv')