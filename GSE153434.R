#-------------------------------------------------------------Figure2: WGCNA of GSE153434----------------------------------------------------------------------------------
#step 0: Data preprocessing
setwd('D:/科研工作站/生信工作站/My project/P003/GSE153434/wgcna')
exp.counts <- read.table(file = 'expression/GSE153434_all.counts.txt',sep = '\t',header = T)
#extract protein coding genes 
exp.counts <- exp.counts[exp.counts$type_of_gene=='protein-coding',c(1,3:22)]
rownames(exp.counts) <- NULL
exp.counts <- exp.counts %>% column_to_rownames('GeneID')
saveRDS(exp.counts,file = 'exp/exp.counts.rds')
#perform log2(1+x) transformation
log21p <- function(x){log2(1+x)}
exp.counts.log <- exp.counts %>% apply(2,log21p) %>% round(2)
#remove stable genes
mad <- apply(exp.counts.log,1,mad)
mad <- order(mad,decreasing = T)
#extract top 5000 genes
exp.counts.log <- exp.counts.log[mad[1:5000],] %>% t() %>% as.data.frame()
#determine whether the genes and samples are qualified
gsg = goodSamplesGenes(exp.counts.log)
gsg$goodSamples %>% table()
gsg$goodGenes %>% table()
#print names and gat rid of bad genes or samples, all samples are good while 717 genes are bad
if (!gsg$allOK){
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(exp.counts.log)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(exp.counts.log)[!gsg$goodSamples], collapse = ", ")));
  exp.counts.log = exp.counts.log[gsg$goodSamples, gsg$goodGenes]
}
#cluster the samples
trait <- data.frame(Patient=c(rep(0,10),rep(1,10)),Normal=c(rep(1,10),rep(0,10)))
rownames(trait) <- rownames(exp.counts.log)
sampleTree <- hclust(dist(exp.counts.log), method = "average")
#transform the trait information into color, white represents low, red represents high and gray represents deficiency
traitColors <-  numbers2colors(trait,signed = F)
#draw samples clustering and group heatmap (Figure2 A)
sizeGrWindow(10,5)
pdf(file = 'md/samples clustering.pdf',width = 8,height = 4)
plotDendroAndColors(sampleTree, traitColors,groupLabels = names(trait), main = "Sample dendrogram and trait heatmap")
abline(h = 120, col = "red")
dev.off()
#remove outliers samples
clust = cutreeStatic(sampleTree,cutHeight = 120, minSize = 2)
exp.counts.log = exp.counts.log[clust%in%c(1,2),]
saveRDS(exp.counts.log,file = 'exp/exp.counts.log.rds')
#step 1:Constructing network step by step
library(WGCNA)
enableWGCNAThreads()
#select soft threshold (Figure2 B)
powers = c(c(1:10), seq(from = 12, to=20, by=2))
RpowerTable=pickSoftThreshold(exp.counts.log, powerVector=powers)[[2]]
write.table(RpowerTable,file = 'soft threshold/soft threshold table.txt',row.names = F,quote = F,sep = '\t')
#plot
pdf(file = 'soft threshold/soft threshold.pdf',width = 8,height = 5)
sizeGrWindow(9,5)
par(mfrow=c(1,2))
plot(RpowerTable[,1], -sign(RpowerTable[,3])*RpowerTable[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",main = paste("Scale independence"))
text(RpowerTable[,1], -sign(RpowerTable[,3])*RpowerTable[,2], labels=powers,cex=0.7,col="red")
abline(h=0.8,col="red")
plot(RpowerTable[,1], RpowerTable[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(RpowerTable[,1], RpowerTable[,5], labels=powers, cex=0.7,col="red")
dev.off()
#The adjacency matrix, topological matrix and dissimilarity matrix are calculated
adjacency <-  adjacency(exp.counts.log,type = 'signed',power= 8)
TOM <-  TOMsimilarity(adjacency)
dissTOM <-  1- TOM
saveRDS(dissTOM,file = 'gene clustering/dissTOM.rds')
#perform hierarchical clustering to genes with the dissimilarity matrix used as distance matrix
library(flashClust)
geneTree <-  flashClust(as.dist(dissTOM),method="average")
#gene modules was built initially and the module number of each gene was returned
dynamicMods <-  cutreeDynamic(dendro = geneTree,distM = dissTOM,deepSplit = 2,pamRespectsDendro = F,minClusterSize = 30)
#transform module numbers into colors
dynamicColors <-  labels2colors(dynamicMods)
#count modules 
unique(dynamicColors) %>% length()
module.gene.counts <- table(dynamicColors) %>% as.data.frame()
write.table(module.gene.counts,file = 'gene clustering/module.gene.counts.txt',sep = '\t',row.names = F,quote = F)
#visualize gene clustering and modules
pdf(file = 'gene clustering/plotDendroAndColors.pdf',width = 8,height = 5)
plotDendroAndColors(dendro = geneTree, colors = dynamicColors, groupLabels = "Dynamic Tree Cut",dendroLabels = F,hang = 0.03,addGuide = T,main = "Gene dendrogram and module colors")
dev.off()
#cluster modules, calculate eigengenes (first principal component) of each modlue
MEList = moduleEigengenes(exp.counts.log,colors=dynamicColors,excludeGrey = T)
MEs = MEList$eigengenes
#calculate inverse correlation
MEDiss1 = 1-cor(MEs)
#module clustering before merging
METree1 = flashClust(as.dist(MEDiss1),method="average")
#merge modules with similarity more than 0.75
merge <-  mergeCloseModules(exp.counts.log, dynamicColors,cutHeight =0.25,verbose=3)
mergedColors <-  merge$colors
module.merge.counts <- table(mergedColors) %>% as.data.frame()
mergedMEs <-  merge$newMEs
#nmber the merged module colors
moduleColors = mergedColors
colorOrder = c("grey",standardColors(50))
moduleLabels = match(moduleColors,colorOrder)-1
#remove gray module, cluster modules after merging
MEs = mergedMEs[,-ncol(MEs)]
MEDiss2 = 1-cor(MEs[,-ncol(MEs)])
METree2 = flashClust(as.dist(MEDiss2),method="average")
#module clustering after merging compared to before
par(mfrow=c(1,2))
plot(METree1,xlab="",sub="",main="Clustering of ME before combined")
abline(h=2.25,col="red")
plot(METree2,xlab="",sub="",main="Clustering of ME after combined")
dev.off()
#gene clustering after module merging compared to before (Figure2 C)
pdf(file = 'gene clustering/plotDendroAndColors(merge).pdf',width = 8,height = 5)
plotDendroAndColors(dendro = geneTree,colors = cbind(dynamicColors,mergedColors),groupLabels = c("Dynamic Tree Cut","Merged Dynamics"),dendroLabels = F, hang = 0.03,addGuide=T,guideHang=0.05,main="Gene Dendrogram and module colors")
dev.off()
#save data
save(MEs,moduleColors,geneTree,file="gene clustering/geneTree.rdata")
#step 2:Module mining
#encode trait information of samples
trait <- data.frame(Patient=c(rep(0,8),rep(1,10)),Normal=c(rep(1,8),rep(0,10)))
#calculate the correlation between modules and traits 
MEs=orderMEs(MEs)
moduleTraitCor <- cor(MEs, trait, use="p")
write.table(moduleTraitCor,file="module mining/module-trait-corr.txt",sep="\t",quote=F)
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor,18)
write.table(moduleTraitPvalue,file="module mining/module-trait-pvalue.txt",sep="\t",quote=F)
#visualization
textMatrix=paste(signif(moduleTraitCor,2),"\n(",signif(moduleTraitPvalue,1),")",sep="")
dim(textMatrix)=dim(moduleTraitCor)
pdf(file = 'module-trait-relationship.pdf',width = 3.5,height = 8)
labeledHeatmap(Matrix=moduleTraitCor,xLabels=colnames(trait),yLabels=names(MEs),ySymbols=names(MEs),
               colorLabels=F,colors=blueWhiteRed(50),textMatrix=textMatrix,setStdMargins=F,cex.text=1,
               cex.lab=1,zlim=c(-1,1),main=paste("Module-trait relationships"))
dev.off()
#correlation between a single module and a certain trait
#extract the name of color of a certain module
modNames = substring(names(MEs), 3)
#calculate the correlation and p value between genes and module
geneModuleMembership = as.data.frame(cor(exp.counts.log, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership),18))
#rename the results
names(geneModuleMembership) <-  paste("MM", modNames, sep="")
names(MMPvalue) <-  paste("p.MM", modNames, sep="")
#calculate the correlation and p value between genes and the objective trait
geneTraitSignificance = as.data.frame(cor(exp.counts.log,trait, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance),18))
#rename the results
names(geneTraitSignificance) <-  paste("GS.", names(trait), sep="")
names(GSPvalue) <-  paste("p.GS.", names(trait)[1], sep="")
save(geneModuleMembership,MMPvalue,geneTraitSignificance,GSPvalue,file = 'MM-GS.rdata')
#set the name of module to be analysed
module = "pink"
#extract genes in the objective module
column = match(module, modNames)
moduleGenes <-  moduleColors==module
#vasualize correlation of the pink module and patient group
pdf(file = 'GS-MM-pink.pdf',width = 5,height = 5)
verboseScatterplot(geneModuleMembership[moduleGenes, column],geneTraitSignificance[moduleGenes,1],
                   xlab = paste("Module Membership"),ylab = "Gene significance",main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
abline(h=0.6,v=0.6)
dev.off()
#select significient genes
mm <- geneModuleMembership[moduleGenes,]
gs <- geneTraitSignificance[moduleGenes,]
gs_mm <- cbind(mm,gs) %>% rownames_to_column('symbol')
writexl::write_xlsx(gs_mm,path = 'pink_gs_mm.xlsx')
pinktop <- pink_gs_mm %>% dplyr::filter(MMpink>0.6,GS.Patient>0.6) %>% .$symbol
write.table(pinktop,file = 'pinktop.txt',sep = '\t',row.names = F,quote = F)
#--------------------------------------------------Figure3 A: volcano plot of DEGs of GSE153434-----------------------------------------------------------------------------
setwd('D:/科研工作站/生信工作站/My project/P003/GSE153434/deg')
deg <- readRDS(file = 'deg.rds')
ggplot(deg,aes(x=logFC,y=-log10(padj),color=change))+
  geom_point()+
  scale_color_manual(values = c('red','blue','gray'))+
  xlim(-10,10)+
  theme_classic()+
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  geom_vline(xintercept = c(-1,1),
             lty=2,
             lwd=0.5,
             col="black")+
  geom_hline(yintercept = -log10(0.05),
             lty=2,
             lwd=0.5,
             col="black")+
  xlab('log2(Fold Change)')+
  ylab('-log10(adj. P value)')+
  labs(title = 'GSE153434')
ggsave(filename = 'volcanoplot_deg.pdf',width = 6,height = 5)
#------------------------------Figure3 B: Venn diagram between DEGs and top genes in pink module of GSE153434----------------------------------------------------
setwd('D:/科研工作站/生信工作站/My project/P003/GSE153434/gene intersection')
deg <- readRDS(file = 'deg.rds')
pinktop <- read.table(file = '../wgcna/module mining/pinktop.txt',header = T)$x
up <- deg$symbol[deg$change=='Up']
geneint <- intersect(up,pinktop) %>% as.data.frame()
write_tsv(geneint,file = 'geneint.tsv')
#VennDiagram
library(VennDiagram)
library(cowplot)
library(ggplotify)
vennplot <- venn.diagram(list('Up-regulated genes'=up,
                              'Pink module top genes'= pinktop),
                         filename = NULL,
                         cex=1.5,
                         col="transparent",
                         fill=c('#AFCEBE','#DD978F'),
                         cat.cex=c(1,1),
                         cat.pos= c(0,0),
                         cat.dist=c(0.05,0.05),
                         scaled = F)
vennplot <- as.ggplot(as_grob(vennplot))
plot(vennplot)
ggsave(filename = 'intersected_genes_up_n_pink.pdf',width = 4,height = 4)
#-----------------------------------------Figure3 E: Boxplot of the hub genes in GSE153434-----------------------------------------------------------------------
#extract deg results of the hub genes
hub <- c('HIF1A','HGF','ITGA5','HMOX1','ITGB3')
deg.hub <- deg[deg$symbol%in%hub & deg$change%in%c('Up','Down'),] %>% arrange(change)
table(deg.hub$change)
writexl::write_xlsx(deg.hub,path = 'deg.hub.xlsx')
##calculate TPM
#load gff3 file
library(GenomicFeatures)
txdb <- makeTxDbFromGFF(file = 'D:/科研工作站/生信工作站/Data base/Genome/Homo_sapiens.GRCh38.112.chr.gff3.gz', format = "gff3")
#get gene length infomation
exons_gene <- exonsBy(txdb, by = "gene")
exons_gene_len <- lapply(exons_gene,function(x){sum(width(reduce(x)))})
exons_gene_len <- as.matrix(t(exons_gene_len))
exons_gene_len <- exons_gene_len[,-1] %>% t() %>% as.data.frame() %>% rownames_to_column('ensid')
colnames(exons_gene_len) <- c('ENSEMBL','exon_len')
write.csv(exons_gene_len,"exons_gene_len.csv", row.names = T)
#integrate expresson matrix
library(clusterProfiler)
library(org.Hs.eg.db)
e2s <- bitr(exons_gene_len$ENSEMBL,fromType = 'ENSEMBL',toType = 'SYMBOL',OrgDb = org.Hs.eg.db)
e2s  <-  e2s[!duplicated(e2s$SYMBOL),]
exons_gene_len <- inner_join(exons_gene_len,e2s,by='ENSEMBL')
exp.counts <- readRDS(file = '../wgcna/expression matrix/exp.counts.rds')
exp.counts <- exp.counts %>% rownames_to_column('SYMBOL')
exp.counts <-  inner_join(exons_gene_len,exp.counts,by='SYMBOL')
exp.counts <- exp.counts[,-1] %>% column_to_rownames('SYMBOL')
#calculate TPM
countToTPM <- function(counts,exon_len){
  A <- counts/exon_len
  sumA <- sum(A)
  exp(log10(A)+log10(1e6)-log10(sumA))
}
exp.tpm <- apply(exp.counts[,2:21],MARGIN = 2,countToTPM,exon_len=exp.counts[,1]) %>% as.data.frame()
saveRDS(exp.tpm,file = 'exp.tpm.rds')
#boxplot
exp.hub <- exp.tpm[hub,] %>% rownames_to_column('Symbol')
exp.hub <- pivot_longer(exp.hub,cols = 2:21,names_to = 'Sample',values_to = 'TPM')
md <- readRDS(file = '../wgcna/meta data/md.rds')
exp.hub$Group <- rep(md$group,5)
ggplot(data = exp.hub,aes(x=Symbol,y=TPM,fill=Group))+geom_boxplot()+theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),legend.position = 'top')+xlab('')+ylab('TPM')+ggtitle('GSE153434')
ggsave(filename = 'boxplot_hubgene.pdf',width = 6,height = 3.5)
#-----------------------------------------Figure7 A: Vocano plot of correlated genes of HIF1A-------------------------------------------------
setwd('D:/科研工作站/生信工作站/My project/P003/GSE153434/corr/HIF1A')
#calculate correlation coefficient
exp.tpm <- readRDS("D:/科研工作站/生信工作站/My project/P003/GSE153434/deg/exp.tpm.rds")
data <- t(exp.tpm)
y <-as.numeric(data[,"HIF1A"])
colnames <-colnames(data)
cor_data_df <-data.frame(colnames)
for(i in 1:length(colnames)){
  print(i)
  test <-cor.test(as.numeric(data[,i]),y,type="spearman")
  cor_data_df[i,2] <-test$estimate
  cor_data_df[i,3] <-test$p.value
}
names(cor_data_df) <-c("Symbol","correlation","pvalue")
cor_data_df <- cor_data_df %>% na.omit()
cor_data_df$direction <- if_else(abs(cor_data_df$correlation)>0.6&cor_data_df$pvalue<0.05,if_else(cor_data_df$correlation>0.6,'Positive','Negative'),'Uncorrelated')
cor_data_df$direction %>% table()
saveRDS(cor_data_df,file = 'cor_data_df.rds')
#-----------------------------------------Figure7 B: Barplot plot of GO enrichment of HIF1A-------------------------------------------------
#GO enrichment for correlated genes
library(clusterProfiler)
library(org.Hs.eg.db)
pos <- cor_data_df$Symbol[cor_data_df$direction=='Positive']
neg <- cor_data_df$Symbol[cor_data_df$direction=='Negative']
go.pos <- enrichGO(pos,OrgDb = org.Hs.eg.db,ont = 'BP',keyType = 'SYMBOL')@result
go.neg<- enrichGO(neg,OrgDb = org.Hs.eg.db,ont = 'BP',keyType = 'SYMBOL')@result
go <- list(pos=go.pos,neg=go.neg)
writexl::write_xlsx(go,path = 'go.xlsx')
#-----------------------------------------Figure8 A: Vocano plot of correlated genes of HGF-------------------------------------------------
setwd('D:/科研工作站/生信工作站/My project/P003/GSE153434/corr/HGF')
#calculate correlation coefficient
exp.tpm <- readRDS("D:/科研工作站/生信工作站/My project/P003/GSE153434/deg/exp.tpm.rds")
data <- t(exp.tpm)
y <-as.numeric(data[,"HGF"])
colnames <-colnames(data)
cor_data_df <-data.frame(colnames)
for(i in 1:length(colnames)){
  print(i)
  test <-cor.test(as.numeric(data[,i]),y,type="spearman")
  cor_data_df[i,2] <-test$estimate
  cor_data_df[i,3] <-test$p.value
}
names(cor_data_df) <-c("Symbol","correlation","pvalue")
cor_data_df <- cor_data_df %>% na.omit()
cor_data_df$direction <- if_else(abs(cor_data_df$correlation)>0.6&cor_data_df$pvalue<0.05,if_else(cor_data_df$correlation>0.6,'Positive','Negative'),'Uncorrelated')
cor_data_df$direction %>% table()
saveRDS(cor_data_df,file = 'cor_data_df.rds')
#-----------------------------------------Figure8 B: Barplot plot of GO enrichment of HGF----------------------------------------------------
#GO enrichment for correlated genes
library(clusterProfiler)
library(org.Hs.eg.db)
pos <- cor_data_df$Symbol[cor_data_df$direction=='Positive']
neg <- cor_data_df$Symbol[cor_data_df$direction=='Negative']
go.pos <- enrichGO(pos,OrgDb = org.Hs.eg.db,ont = 'BP',keyType = 'SYMBOL')@result
go.neg<- enrichGO(neg,OrgDb = org.Hs.eg.db,ont = 'BP',keyType = 'SYMBOL')@result
go <- list(pos=go.pos,neg=go.neg)
writexl::write_xlsx(go,path = 'go.xlsx')