##step 0:数据预处理--------------
##整理表达矩阵
library(WGCNA)
raw.counts <- read.table(file = 'GSE153434\\series_matrix\\GSE153434_all.counts.txt',sep = '\t',header = T)
table(raw.counts$type_of_gene)
#提取PCG
counts <- raw.counts[raw.counts$type_of_gene=='protein-coding',c(1,3:22)]#得到5310个
rownames(counts) <- NULL
counts <- counts %>% column_to_rownames('GeneID')
write.table(counts,file = 'GSE153434\\series_matrix\\GSE153434_counts.txt',quote = F,sep = '\t',row.names = T)
#log2(1+x)转化
log21p <- function(x){log2(1+x)}
counts.log <- counts %>% apply(2,log21p) %>% round(2)
write.table(counts.log,file = 'GSE153434\\series_matrix\\GSE153434_counts_log.txt',quote = F,sep = '\t',row.names = T)
##去除变化小的基因
#计算基因绝对中位差，降序排列
count.log_mad <- apply(counts.log,1,mad)
count.log_mad_sorted <- order(count.log_mad,decreasing = T)
count.log_mum <- count.log_mad_sorted[1:5000]
#提取前5000个基因
count.log_filter <- counts.log[count.log_mum,] %>% t() %>% as.data.frame()
#保存
save(count.log_filter,file = 'GSE153434\\series_matrix\\count.log_filter.rdata')
##去除低表达的样本和基因
datExpr0 <- count.log_filter
gsg = goodSamplesGenes(datExpr0)#判断基因样本是否合格，默认一个基因至少在一半样本中表达
gsg$goodSamples %>% table()
#TRUE 
#20
gsg$goodGenes %>% table()
#FALSE  TRUE 
#717  4283 #20个样本都合格，717个基因不合格
#去除基因和样本，打印不合格的基因和样本名称
if (!gsg$allOK){
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # 提取保留的基因和样本
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
##样本聚类去除异常样本
sampleTree = hclust(dist(datExpr0), method = "average")#计算样本间的距离，类平均法
sizeGrWindow(10,5)
pdf(file = "GSE153434\\plot\\sample_cluster.pdf", width = 10, height = 5);
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)
#根据上图判断，需要截取的高度参数h
abline(h = 120, col = "red")#在120的地方画条线
dev.off()
#去除离群得聚类样本，cutHeight参数要与上述得h参数值一致，
clust = cutreeStatic(sampleTree, cutHeight = 120, minSize = 2)
table(clust)
#0  1  2 
#2 10  8 
#提取保留样本，去除了对照组的4和10，最终剩10个疾病，8个对照
keepSamples = (clust%in%c(1,2))
datExpr = datExpr0[keepSamples,]
#基因样本数
nGenes = ncol(datExpr);nSamples = nrow(datExpr)
#4328，18
save(datExpr,file = 'GSE153434\\series_matrix\\datExpr.rdata')
##step 1:多步法构建网络 --------------
rm(list = ls())
library(WGCNA)
enableWGCNAThreads()#打开多线程
##测试软阈值
powers = c(c(1:10), seq(from = 12, to=20, by=2))#设置测试的范围
RpowerTable=pickSoftThreshold(datExpr, powerVector=powers)[[2]]
write.table(RpowerTable,file = 'D:\\科研工作站\\生信工作站\\sskt5.0\\GSE153434\\wgcna\\preprocess\\power_table.txt',
            row.names = F,quote = F,sep = '\t')
sizeGrWindow(9,5)
par(mfrow=c(1,2))
plot(RpowerTable[,1], -sign(RpowerTable[,3])*RpowerTable[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",main = paste("Scale independence"))
text(RpowerTable[,1], -sign(RpowerTable[,3])*RpowerTable[,2], labels=powers,cex=0.7,col="red")
abline(h=0.8,col="red")
plot(RpowerTable[,1], RpowerTable[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(RpowerTable[,1], RpowerTable[,5], labels=powers, cex=0.7,col="red")
#选择R^2大于0.8后的最小β，所以选择β=8

##计算邻接矩阵、拓扑矩阵、相异度矩阵
adjacency <-  adjacency(datExpr,type = 'signed',power= 8)#有符号网络
TOM <-  TOMsimilarity(adjacency)
dissTOM <-  1- TOM
##构建基因层次聚类树
library(flashClust)
geneTree <-  flashClust(as.dist(dissTOM),method="average")#使用相异度矩阵作为距离矩阵
#构建初步基因模块
dynamicMods <-  cutreeDynamic(dendro = geneTree,
                            distM = dissTOM,
                            deepSplit = 2,#值越大，模块的数量越多，但每个模块的规模越小
                            pamRespectsDendro = F,#不懂
                            minClusterSize = 30)#设定模块最小基因数
# 将数字标签转换为颜色
dynamicColors <-  labels2colors(dynamicMods)
#模块统计
unique(dynamicColors) %>% length()#14个模块
module.stat <- table(dynamicColors) %>% as.data.frame()#各模块所含基因数
write.table(module.stat,file = 'GSE153434\\plot\\module.stat.txt',sep = '\t',row.names = F,quote = F)
#可视化
plotDendroAndColors(dendro = geneTree, 
                    colors = dynamicColors, 
                    groupLabels = "Dynamic Tree Cut",
                    dendroLabels = F,
                    hang = 0.03,
                    addGuide = T,
                    main = "Gene dendrogram and module colors")
dev.off()
##模块聚类与合并
#计算模块特征基因，即每个模块的第一主成分
MEList = moduleEigengenes(datExpr,colors=dynamicColors,excludeGrey = T)
MEs = MEList$eigengenes
MEDiss1 = 1-cor(MEs)#计算倒相关性距离
#合并前模块聚类
METree1 = flashClust(as.dist(MEDiss1),method="average")
##合并相似度大于0.75的模块
merge <-  mergeCloseModules(datExpr, 
                          dynamicColors,
                          cutHeight =0.25, 
                          verbose=3)
mergedColors <-  merge$colors
table(dynamicColors)#合并前的模块颜色
module.merge.stat <- table(mergedColors) %>% as.data.frame()
write.table(module.merge.stat,file = 'GSE153434\\plot\\module.merge.stat.txt',sep = '\t',row.names = F,quote = F)
#合并后的模块颜色，可以看到从14个模块变成了9个模块
mergedMEs <-  merge$newMEs#合并后9个模块的特征向量
#对合并后的模块颜色编号
moduleColors = mergedColors
colorOrder = c("grey",standardColors(50))
moduleLabels = match(moduleColors,colorOrder)-1
##对合并后的模块进行聚类
MEs = mergedMEs[,-ncol(MEs)]
MEDiss2 = 1-cor(MEs[,-ncol(MEs)])
METree2 = flashClust(as.dist(MEDiss2),method="average")
#模块合并前后聚类比较
par(mfrow=c(1,2))
plot(METree1,xlab="",sub="",main="Clustering of ME before combined")
abline(h=2.25,col="red")#相异度为0.25
plot(METree2,xlab="",sub="",main="Clustering of ME after combined")
dev.off()
#模块合并前后基因聚类树
plotDendroAndColors(dendro = geneTree,#剪切树
                    colors = cbind(dynamicColors,mergedColors),#合并前后的模块比较
                    groupLabels = c("Dynamic Tree Cut","Merged Dynamics"),
                    dendroLabels = F, 
                    hang = 0.03,
                    addGuide=T,
                    guideHang=0.05,
                    main="Gene Dendrogram and module colors")
dev.off()
#保存数据
save(MEs, moduleLabels, moduleColors, geneTree, file="GSE153434\\wgcna\\geneTree.rdata")
##step 2:模块挖掘--------------
##样本与分组
#读取表型信息
trait <- c(rep('Normal',8),rep('TAAD',10)) %>% as.factor() %>% as.numeric() %>% as.data.frame()
rownames(trait) <- rownames(datExpr)
colnames(trait)[1] <- 'TAAD'
# 对样本进行聚类
sampleTree <- hclust(dist(datExpr), method = "average")
#将临床信息转换为颜色，白色表示低，红色表示高，灰色表示缺失
traitColors <-  numbers2colors(trait, signed = F)
#样本聚类与分组热图
plotDendroAndColors(sampleTree, 
                    traitColors,
                    groupLabels = names(trait), 
                    main = "Sample dendrogram and trait heatmap")
dev.off()
#模块表型相关性-----
MEs=orderMEs(MEs)
moduleTraitCor <- cor(MEs, trait, use="p")#处理缺失值的方法"pairwise.complete.obs"
write.table(moduleTraitCor,file="GSE153434\\wgcna\\modtrait.cor.txt",sep="\t",quote=F)
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor,18)
write.table(moduleTraitPvalue,file="GSE153434\\wgcna\\modtrait.p.txt",sep="\t",quote=F)
#可视化
textMatrix=paste(signif(moduleTraitCor,2),"\n(",signif(moduleTraitPvalue,1),")",sep="")
dim(textMatrix)=dim(moduleTraitCor)
#基因模块与临床信息相关性图
pdf(file = 'GSE153434\\plot\\Module_trait_relationship.pdf',width = 3.5,height = 8)
labeledHeatmap(Matrix=moduleTraitCor,
               xLabels=colnames(trait),
               yLabels=names(MEs),
               ySymbols=names(MEs),
               colorLabels=F,
               colors=blueWhiteRed(50),
               textMatrix=textMatrix,
               setStdMargins=F,
               cex.text=1,
               cex.lab=1,
               zlim=c(-1,1),
               main=paste("Module-trait relationships"))
dev.off()
##单一模块与某一表型相关性
# 模块对应的颜色
modNames = substring(names(MEs), 3)#提取模块的颜色名
# 计算模块与基因的相关系数、P值
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership),18))
# 对结果进行命名
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

# 计算目标表型与基因的相关系数、P值
geneTraitSignificance = as.data.frame(cor(datExpr,trait$TAAD, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance),18))
# 对结果进行命名
names(geneTraitSignificance) = paste("GS.", names(trait), sep="")
names(GSPvalue) = paste("p.GS.", names(trait), sep="")
save(geneModuleMembership,MMPvalue,geneTraitSignificance,GSPvalue,file = 'GSE153434\\wgcna\\MM_GS.rdata')
# 设置需要分析的模块名称，此处为brown模块
module = "pink"
module = "green"
# 提取brown模块中的基因
column = match(module, modNames)
moduleGenes <-  moduleColors==module
# 可视化brown模块与M分期的相关性分析结果
pdf(file = paste0('GSE153434\\plot\\',module,'_trait_cor.pdf'),width = 5,height = 5)
verboseScatterplot(geneModuleMembership[moduleGenes, column],
                   geneTraitSignificance[moduleGenes, 1],
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for TAAD",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()
#模块枢纽基因----
#一个模块选出一个枢纽基因
HubGenes <- chooseTopHubInEachModule(datExpr,moduleColors)
write.table(HubGenes,file = 'GSE153434\\wgcna\\hubgenes_modules.txt')
#表型相关的枢纽基因
NS = networkScreening(trait$TAAD,MEs,datExpr)
write.table(NS,file = 'GSE153434\\wgcna\\hubgenes_TAAD.txt')
#模块内富集----
library(anRichment)
library(clusterProfiler)
#GO富集
#构建GO背景基因集
GOcollection = buildGOcollection(organism = "human")
geneNames = colnames(datExpr)
# 将基因SYMBOL转换为ENTREZID基因名
geneID = bitr(geneNames,fromType = "SYMBOL", toType = "ENTREZID", 
              OrgDb = "org.Hs.eg.db", drop = FALSE)
write.table(geneID, file = "GSE153434\\wgcna\\geneID_map.txt", sep = "\t", quote = F, row.names = F)
# 进行GO富集分析
GOenr = enrichmentAnalysis(classLabels = moduleColors,
                           identifiers = geneID$ENTREZID,
                           refCollection = GOcollection,
                           useBackground = "given",
                           threshold = 1e-4,
                           thresholdType = "Bonferroni",
                           getOverlapEntrez = TRUE,
                           getOverlapSymbols = TRUE,
                           ignoreLabels = "grey");
#保存结果
tab = GOenr$enrichmentTable
write.table(tab, file = "GSE153434\\wgcna\\GOEnrichmentTable.txt", sep = "\t", quote = F, row.names = F)

# 提取主要结果，并写入文件
screenTab <- tab %>% dplyr::select(1, 3, 4, 6, 7, 8, 13)
screenTab[, c(4, 5, 6)] <- screenTab[, c(4, 5, 6)] %>% apply(2,as.numeric) %>% signif(2)#简洁
# 给结果命名
colnames(screenTab) = c("module", "GOID", "term name", "p-val", "Bonf", "FDR", "size")
rownames(screenTab) = NULL

# 查看结果
head(screenTab)
# 写入文件中
write.table(screenTab, file = "GSE153434\\wgcna\\GOEnrichmentTableScreen.txt", sep = "\t", quote = F, row.names = F)

#KEGG富集分析
# AnRichment没有直接提供KEGG数据的背景集,这里使用MSigDBCollection构建C2通路数据集
KEGGCollection = MSigDBCollection("GSE153434\\msigdb_v7.1.xml",
                                  MSDBVersion = "7.1",
                                  organism = "human",
                                  excludeCategories = c("h","C1","C3","C4","C5","C6","C7")) 
# KEGG分析
KEGGenr = enrichmentAnalysis(classLabels = moduleColors,
                             identifiers = geneID$ENTREZID,
                             refCollection = KEGGCollection,
                             useBackground = "given",
                             threshold = 1e-4,
                             thresholdType = "Bonferroni",
                             getOverlapEntrez = TRUE,
                             getOverlapSymbols = TRUE,
                             ignoreLabels = "grey")
# 提取KEGG结果，并写入文件
tab = KEGGenr$enrichmentTable
names(tab)
write.table(tab, file = "GSE153434\\wgcna\\KEGGEnrichmentTable.txt", sep = "\t", quote = F, row.names = F)

# 提取主要结果并写入文件
keepCols = c(1, 3, 4, 6, 7, 8, 13)
screenTab = tab[, keepCols]
# 取两位有效数字
numCols = c(4, 5, 6)
screenTab[, numCols] = signif(apply(screenTab[, numCols], 2, as.numeric), 2)

# 对结果表格进行重命名
colnames(screenTab) = c("module", "ID", "term name", "p-val", "Bonf", "FDR", "size")
rownames(screenTab) = NULL
# 查看结果
head(screenTab)
# 写入文件中
write.table(screenTab, file = "GSE153434\\wgcna\\KEGGEnrichmentTableScreen.txt", sep = "\t", quote = F, row.names = F)
#富集分析可视化
library(openxlsx)
GOEnrichmentTable <- read.xlsx('GSE153434\\wgcna\\enrichment\\GOEnrichmentTable.xlsx')
colnames(GOEnrichmentTable)
module <- 'pink'
table <- GOEnrichmentTable[GOEnrichmentTable$class%in%module,]
pdf(file = paste0('GSE153434\\wgcna\\enrichment\\',module,'_barplot.pdf'))
ggplot(table,aes(-log10(FDR),reorder(dataSetName,-log10(FDR))))+
  geom_bar(stat = 'identity',width = 0.6,fill=module)+theme_classic()+
  theme(legend.position = 'none')+
  ylab('')
dev.off()

##step 3:输出cytoscape可视化-----------
#重新计算TOM,power值设置为前面选择好的
adjacency <-  adjacency(datExpr,type = 'signed',power= 8)
TOM <-  TOMsimilarity(adjacency)
##输出全部网络模块
cyt = exportNetworkToCytoscape(TOM,
                               edgeFile = "GSE153434\\wgcna\\CytoscapeInput-edges-all.txt",
                               nodeFile = "GSE153434\\wgcna\\CytoscapeInput-nodes-all.txt",
                               weighted = T,
                               threshold = 0.1,
                               nodeNames = geneID$SYMBOL,
                               altNodeNames = geneID$ENTREZID,
                               nodeAttr = moduleColors)
# 输出感兴趣网络模块
#提取目标模块的基因名字和TOM子矩阵，保存为网络文件
modules = "green"
modules = "pink"
#提取目标模块的基因
inModule = is.finite(match(moduleColors, modules))#返回第一集合中的每一个元素在第二个集合中的位次
modGenes = geneID[inModule,]
# 选择指定模块的TOM矩阵
modTOM = TOM[inModule, inModule]

# 输出为Cytoscape软件可识别格式
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("GSE153434\\wgcna\\CytoscapeInput-edges-green", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("GSE153434\\wgcna\\CytoscapeInput-nodes-green", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = T,
                               threshold = 0.2,
                               nodeNames = modGenes$SYMBOL,
                               altNodeNames = modGenes$ENTREZID,
                               nodeAttr = moduleColors[inModule])


##step 4:差异基因分析-----
#GSE98770----
#提取表达矩阵，提取表型信息
library(GEOquery)
GSE98770_sm <- getGEO(filename = "DEG\\Series Matrix\\GSE98770-GPL14550_series_matrix.txt.gz")
GSE98770_pd0 <- pData(GSE98770_sm)
GSE98770_pd <- data.frame(sample_id=GSE98770_pd0$geo_accession,
                          sample_name=c(paste0("ATAAD","_",1:6),paste0("Normal","_",1:5)),
                          group=factor(c(rep("ATAAD",6),rep("Normal",5)),levels = c("Normal","ATAAD")),
                          gender=str_sub(GSE98770_pd0$characteristics_ch1.1,8,-1),
                          age=str_sub(GSE98770_pd0$characteristics_ch1.2,5,-2))

#读取原始数据
library(limma)
GSE98770_raw <- read.maimages(files = list.files("DEG\\RAW Data\\txt"),
                              source = "agilent",
                              path = "GSE98770/RAW data",
                              names = GSE98770_targets$sample_id,
                              other.columns = "gIsWellAboveBG",
                              green.only = T)
GSE98770_expr_raw <- GSE98770_raw$E %>% as.data.frame() %>% set_names(GSE98770_pd$sample_id)
#标准化前质量评估
boxplot(GSE98770_expr_raw,
        las=2,
        outline=F,
        col=c(rep("#FF6666",6),rep("#9999FF",5)),
        ylab="Expression Value",
        xlab="Sample List",
        names=F,
        main="GSE98770 before normalization")
legend("topright",
       fill = c("#FF6666","#9999FF"),
       legend = c("ATAAD","Control"))

#背景校正和标准化
library(limma)
GSE98770_bgc <- limma::backgroundCorrect(RG=GSE98770_RAW,
                                         method = "normexp",
                                         offset = 50,#补偿值
                                         normexp.method = "mle")
GSE98770_norm <- limma::normalize=BetweenArrays(GSE98770_bgc,
                                                method = "quantile")
GSE98770_expr_norm <- GSE98770_norm$E
rownames(GSE98770_expr_norm) <- GSE98770_norm$genes$ProbeName
#校正后质量评估
boxplot(GSE98770_expr_norm,
        las=2,
        outline=F,
        col=c(rep("#FF6666",6),rep("#9999FF",5)),
        ylab="Expression Value",
        xlab="Sample List",
        names=F,
        main="GSE98770 after normalization")
legend("topright",
       fill = c("#FF6666","#9999FF"),
       legend = c("ATAAD","Control"))
#arrayQM质量评估
library(arrayQualityMetrics)
rownames(GSE98770_pd)<- colnames(GSE98770_norm$E) <- GSE98770_pd$sample_id
GSE98770_est <- ExpressionSet(assayData =GSE98770_norm$E,
                              phenoData = AnnotatedDataFrame(data = GSE98770_pd))
arrayQualityMetrics(GSE98770_est,
                    outdir = "GSE98770/arrayQM",
                    force = T,
                    intgroup = "group",
                    do.logtransform = T)
#芯片注释
source("annotated_expr.R")
GSE98770_p2s <- read_tsv("DEG\\Annotation\\028004_D_GeneList_20210927.txt") %>% 
  dplyr::select(1,3) %>% set_names("probe_id","symbol")

GSE98770_expr_anno <- annotate_expr(expr = GSE98770_expr_norm,p2s = GSE98770_p2s) %>%
  t() %>% as.data.frame() %>% 
  rownames_to_column("sample_id") %>%
  inner_join(GSE98770_pd,.,by="sample_id") %>% 
  arrange(group) %>% 
  column_to_rownames("sample_id") %>% 
  dplyr::select(-1,-2,-3) %>% t() %>% as.data.frame()

#差异分析
expr <- GSE98770_expr_anno
pd <- GSE98770_pd
#用sva探索批次效应并计算替代变量
library(sva)
#建立两个设计矩阵，一个含分组，一个只有截距
mod <- model.matrix(~group,data = pd)
mod0 <- model.matrix(~1,data = pd)
#用num.sv估计批次效应的个数
n.sv <- num.sv(dat=expr,mod,method = "leek")
#用sva计算每个样本4个批次效应的替代变量，加到mod中去
svobj <- sva(as.matrix(expr ),mod=mod,mod0 = mod0,n.sv = n.sv)
modSv <- cbind(mod,svobj$sv)
#纳入未知的批次效应
library(limma)
design_sva <-modSv
aw_sva <- arrayWeights(expr,design = design_sva)
fit_sva <- lmFit(expr, design_sva,weights = aw_sva)%>% eBayes(trend = T)
res_DEG_sva<- topTable(fit_sva,coef =2,n = Inf)
res_DEG_sva_sig<- topTable(fit_sva,coef =2,n = Inf,p=0.05,lfc=0)
GSE98770_DEG <- res_DEG_sva
#标记上下调
lfc <- 1;p <- 0.05
up <- (GSE98770_DEG$logFC>lfc)&(GSE98770_DEG$adj.P.Val<p)
dw <- (-GSE98770_DEG$logFC>lfc)&(GSE98770_DEG$adj.P.Val<p)
change <- ifelse(up,'up',ifelse(dw,'down','stable'))
table(change)
GSE98770_DEG$change <- change
###保存
save(GSE98770_pd,GSE98770_expr_anno,GSE98770_DEG,file = "GSE98770_DEG.rdata")
save(GSE98770_raw,GSE98770_norm,file = "GSE98770_process.rdata")
#火山图
library(tidyverse)
library(ggplot2)
library(ggrepel)
DEG_all <- GSE98770_DEG%>% rownames_to_column("symbol") %>%
  dplyr::select(symbol,logFC,Pvalue=adj.P.Val) %>% 
  mutate(direction=factor(ifelse(Pvalue<0.05&abs(logFC)>1,
                                 ifelse(logFC>1,"Up-regulate","Down-regulate"),
                                 "No-different"),
                          levels =c("Up-regulate","Down-regulate","No-different")))
ggplot(data = DEG_all,aes(x=logFC,y=-log10(Pvalue),col=direction))+
  geom_point(alpha=0.6)+
  scale_color_manual(values=c("red","blue","#808080"))+
  #geom_text_repel(data = DEG_all %>% filter(Pvalue<0.05,abs(logFC)>4),
  #aes(label=symbol),size=3,segment.color="red",show.legend = F)+
  theme_classic()+
  theme(legend.title = element_blank(),
        legend.position = "top")+
  ylab(expression(-log[10]("P Value")))+
  xlab(expression(log[2]("Fold Changge")))+
  xlim(-4,4)+
  ylim(0,6)+
  geom_vline(xintercept = c(-1,1),
             lty=2,
             lwd=0.6,
             col="orange")+
  geom_hline(yintercept = -log10(0.05),
             lty=2,
             lwd=0.6,
             col="orange")
#GSE153434----
#GSE153434使用GEO2R的结果
GSE153434_GEO2R <- read_tsv(file = 'GSE153434\\dge\\GSE153434.top.table.tsv')
dgeall <- GSE153434_GEO2R[,c(8,6,2)] %>% na.omit()
colnames(dgeall) <- c('symbol','logFC','adj.P.Val')
lfc <- 1;p <- 0.05
up <- (dgeall$logFC>lfc)&(dgeall$adj.P.Val<p)
dw <- (-dgeall$logFC>lfc)&(dgeall$adj.P.Val<p)
dgeall$change <- ifelse(up,'up',ifelse(dw,'down','stable'))
dgeall$change <- factor(dgeall$change,levels = c('up','down','stable'))
table(dgeall$change)
type <- read.table('GSE153434\\series_matrix\\GSE153434_all.counts.txt',sep = '\t') %>%.[2:nrow(type),1:2]
dgeall <- dgeall %>% inner_join(type,by='symbol')
colnames(type) <- c('symbol','type')
table(dgeall$type)
dgeall <- dgeall[dgeall$type%in%'protein-coding',]
write.table(dgeall,file = 'GSE153434\\dge\\GSE153434.geo2r.txt',row.names = F,quote = F,sep = '\t')
dgeall <- read.table(file = 'dge\\GSE153434.geo2r.txt',header = T)
dgeall$change <- factor(dgeall$change,levels = c('up','down','stable')) 
#GSE98770
dgeall <- GSE98770_DEG
#差异基因可视化-火山图
pdf(file = 'dge\\GSE98770_valcano.pdf',width = 8,height = 4)
pdf(file = 'dge\\GSE153434_valcano.pdf',width = 8,height = 4)
ggplot(data = dgeall,aes(x=logFC,y=-log10(adj.P.Val),col=change))+
  geom_point(alpha=0.6)+
  scale_color_manual(values=c("red","blue","gray"))+
  theme_classic()+
  theme(legend.title = element_blank(),
        legend.position = "top")+
  ylab(expression(-log[10]("P Value")))+
  xlab(expression(log[2]("Fold Changge")))+
  xlim(-5,5)+
  ylim(0,6)+
  geom_vline(xintercept = c(-1,1),
             lty=2,
             lwd=0.6,
             col="orange")+
  geom_hline(yintercept = -log10(0.05),
             lty=2,
             lwd=0.6,
             col="orange")
dev.off()
#差异基因可视化-热图
library(pheatmap)
library(RColorBrewer)
counts.log <- read.table('GSE153434\\series_matrix\\GSE153434_counts_log.txt')
counts.log <- counts.log[,-c(2,5)]
genes <- dgelist[top[1:1000],] %>% rownames()
heatmap_data <- counts[rownames(counts.log)%in%genes,]
annotation_col <- data.frame(sample=colnames(counts),
                             group=c(rep('Normal',8),rep('Patient',10))) %>% column_to_rownames('sample')
pdf(file = 'GSE153434\\expression\\heatmap.pdf',width = 6,height = 8)
pheatmap(heatmap_data,
         cluster_rows = F,
         cluster_cols = F,
         scale = 'row',
         annotation_col = annotation_col,
         annotation_colors = list(group=c(Normal='blue',Patient='red')),
         show_colnames = F,
         show_rownames = F)
dev.off()
#差异基因取交集
GSE153434_DEG <- read.table('GSE153434\\dge\\GSE153434.geo2r.txt',header = T)
GSE153434_DEG <- GSE153434_DGE[abs(GSE153434_DGE$logFC)>1&GSE153434_DGE$adj.P.Val<0.05,]
GSE98770_DEG <- GSE98770_DEG[abs(GSE98770_DEG$logFC)>1&GSE98770_DEG$adj.P.Val<0.05,]
GSE98770_dge <- rownames(GSE98770_DGE)
GSE153434_dge <- GSE153434_DGE$symbol
int <- intersect(GSE153434_dge,GSE98770_dge)
GSE98770_single <- setdiff(GSE98770_dge,int)
GSE153434_single <- setdiff(GSE153434_dge,int)
#可视化
library(VennDiagram)
library(cowplot)
library(ggplotify)
vennplot <- venn.diagram(list("GSE153434"=GSE153434_dge,
                              "GSE98770"=GSE98770_dge),
                           filename = NULL,
                           cex=1.5,
                           col="transparent",
                           fill=c('skyblue','gray'),
                           cat.cex=c(1,1),
                           cat.pos= c(0,0),
                           cat.dist=c(0.2,0.2),
                           scaled = T)
class(vennplot)
vennplot <- as.ggplot(as_grob(vennplot))
plot(vennplot)
plot(vennplot)
#差异基因功能富集
library(clusterProfiler)
library(org.Hs.eg.db)
GSE153434_dge <- read.table('GSE153434\\dge\\GSE153434.geo2r.txt',header = T)
gene_up <- GSE153434_dge$symbol[GSE153434_dge$change=='up']
gene_dw <- GSE153434_dge$symbol[GSE153434_dge$change=='down']
ego_up <- enrichGO(gene_up,
                   OrgDb = org.Hs.eg.db,
                   keyType = 'SYMBOL',
                   ont = 'ALL',
                   readable = T)
ego_up.df <- ego_up@result
write.table(ego_up.df,file = 'GSE153434\\ego\\ego_up.txt',sep = '\t',quote = F,row.names = F)
ego_dw <- enrichGO(gene_dw,
                   OrgDb = org.Hs.eg.db,
                   keyType = 'SYMBOL',
                   ont = 'ALL',
                   readable = T)
ego_dw.df <- ego_dw@result
write.table(ego_dw.df,file = 'GSE153434\\ego\\ego_dw.txt',sep = '\t',quote = F,row.names = F)
ego_int <- enrichGO(int,
                   OrgDb = org.Hs.eg.db,
                   keyType = 'SYMBOL',
                   ont = 'ALL',
                   readable = T)
ego_int.df <- ego_int@result
write.table(ego_int.df,file = 'GSE153434\\ego\\ego_int.txt',sep = '\t',quote = F,row.names = F)

ego_gse98770_sin <- enrichGO(GSE98770_single,
                    OrgDb = org.Hs.eg.db,
                    keyType = 'SYMBOL',
                    ont = 'ALL',
                    readable = T)
ego_gse98770_sin.df <- ego_gse98770_sin@result
write.table(ego_gse98770_sin.df,file = 'GSE153434\\ego\\ego_gse98770_sin.txt',sep = '\t',quote = F,row.names = F)

ego_gse153434_sin <- enrichGO(GSE153434_single,
                             OrgDb = org.Hs.eg.db,
                             keyType = 'SYMBOL',
                             ont = 'ALL',
                             readable = T)
ego_gse153434_sin.df <- ego_gse153434_sin@result
write.table(ego_gse153434_sin.df,file = 'GSE153434\\ego\\ego_gse153434_sin.txt',sep = '\t',quote = F,row.names = F)
#可视化
pdf(file = 'GSE153434\\ego\\ego_up.pdf',width = 6.5,height = 3)
ggplot(ego_up.df[1:10,],aes(p.adjust,reorder(Description,log10(p.adjust))))+
  geom_bar(stat = 'identity',width = 0.7,fill='red')+theme_classic()
dev.off()
#两个数据GO上调取交集
library(VennDiagram)
library(RColorBrewer)
#PPI网络基因表达----
##pink_module所有基因34个,7个无连接，还剩27个
genes <- 'ADM
ADORA2B
ARG2
GINS1
GINS3
GINS4
HGF
HIF1A
HMGA1
HMOX1
IFNGR1
IL18R1
IL1R1
IL1RAP
IL1RL1
IL1RL2
ITGA2
ITGA5
ITGB3
ITGB6
KLKB1
LBP
S100A8
SERPINA1
SERPINA3
SLC11A1
SNAI2'
genes <- str_split(genes,pattern = '\\s+',simplify = T) %>% as.vector() %>% unique()
GSE153434_DEG <- read.table('GSE153434\\dge\\GSE153434.geo2r.txt',header = T) %>% column_to_rownames('symbol')
dgeall <- GSE153434_DEG
ppi_target <- dgeall[rownames(dgeall)%in%genes,] %>% rownames_to_column('symbol')
write.table(ppi_target,file = 'GSE153434\\ppi\\node.table.txt',row.names = F,sep = '\t',quote = F)
#转录调控网络基因表达----
#靶基因
genes <-'ITGA5
KRT16
ITGB3
HMGA1
HGF
IFNGR1
ITGA2
ALOX15B
HMOX1
ADM
HIF1A
SERPINA3'
genes <- str_split(genes,pattern = '\\s+',simplify = T) %>% as.vector()
#差异表达
tf_target <- dgeall[rownames(dgeall)%in%genes,] %>% rownames_to_column('symbol')#没有ALOX15B，过滤掉了
write.table(tf_target,file = 'GSE153434\\trust\\tf_target.txt',row.names = F,sep = '\t',quote = F)
#表达热图
library(pheatmap)
library(RColorBrewer)
counts <- read.table(file = 'series_matrix\\GSE153434_counts_log.txt')
counts <- counts[,-c(2,5)]
colnames(counts) <- c(paste('Control',c(1,2,3,5,6,7,8,9)),paste('Patient',c(10,1,2,3,4,5,6,7,8,9)))
heatmap_data <- counts[rownames(counts)%in%genes,]
annotation_col <- data.frame(sample=colnames(counts),
                             group=c(rep('Control',8),rep('Patient',10))) %>% column_to_rownames('sample')
pdf(file = 'expression\\heatmap.pdf',width = 8,height = 4)
pheatmap(heatmap_data,
         cluster_rows = F,
         cluster_cols = F,
         scale = 'row',
         annotation_col = annotation_col,
         annotation_colors = list(group=c(Control='blue',Patient='red')),
         show_colnames = T,
         angle_col = 45)
dev.off()
#5个转录因子在GSE98770中的差异表达情况
#转录因子
tf <- 'SP1
TFAP2A
STAT3
RELA
NFKB1'
tf <- str_split(tf,pattern = '\\s+',simplify = T)
tf_dge <- GSE98770_DEG[rownames(GSE98770_DEG)%in%tf,]
pdf(file = 'expression\\tf_expression.pdf',width = 5,height = 4)
ggplot(tf_dge,aes(reorder(rownames(tf_dge),logFC),logFC,fill=-log10(adj.P.Val)))+
  geom_bar(stat = 'identity',width = 0.5)+
  scale_fill_gradient(low='blue',high = 'red')+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))+
  ylab('')
dev.off()
##step 5:ROC分析--------------
library(pROC)
genes <- c('RELA','STAT3','SERPINA3','HIF1A','HMOX1','HGF')
roc.data <- GSE98770_expr_anno[rownames(GSE98770_expr_anno)%in%genes,] %>% t() %>% as.data.frame() 
roc.data$group <- GSE98770_pd$group
roc_list <- list()
for(i in 1:6){
  roc <- roc.data[,c(genes[i],'group')]
  colnames(roc) <- c('gene','group')
  roc_list[[i]] <- roc
}
names(roc_list) <- genes
roc.ob <- list()
for(i in 1:6){
  roc.ob[[i]] <- roc(roc_list[[i]]$group,roc_list[[i]]$gene)
}
for(i in 1:6){
  pdf(file = paste0('GSE153434\\roc\\',genes[i],'_roc.pdf'))
  roc.plot <- plot(roc.ob[[i]],print.auc=T)
  print(roc.plot)
  dev.off()
}
save(roc.data,roc_list,roc.ob,roc.plot,file = 'GSE153434\\roc\\roc_analysis.rdata')