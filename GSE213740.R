#------------------------------------------Figure5: Clustering and dimention reduction-----------------------------------------------------------
#step 1:create a seurat object
setwd('D:/科研工作站/生信工作站/My project/P003/GSE213740')
library(Seurat)
library(BPCells)
samples.dir <- list.files(path = 'D:/科研工作站/生信工作站/My project/P002/GSE213740/exp.mat/matrix',full.names = T)[-1]
#BPC transformation 
data.list <- c()
samples.name <- c(paste0('AD',1:5),paste0('Ctrl',1:3))
for(i in 1:length(samples.name)){
  CreateSeuratObject(counts=Read10X(data.dir = samples.dir[i]))
  x <- RenameCells(x,add.cell.id=samples.name[i])  
  write_matrix_dir(mat=x[['RNA']]$counts,dir = paste0('BPCell',samples.name[i]))
  data.list[[i]] <- open_matrix_dir(dir = paste0('BPCell',samples.name[i]))
}
names(data.list) <- samples.name
saveRDS(data.list,file = 'data.list.rds')
#初始细胞数
for(i in 1:length(data.list))print(paste(names(data.list[i]),ncol(data.list[[i]])))
data.merge <- CreateSeuratObject(counts = data.list,min.cells = 3,min.features = 200)
md <- data.merge@meta.data
md$group <- sapply(strsplit(as.character(md$orig.ident),split = '-'),'[',1)
data.merge@meta.data <- md
data.merge <- AddMetaData(data.merge,metadata = PercentageFeatureSet(data.merge, pattern = "^MT-"),col.name = 'percent.mt')
table(data.merge$orig.ident,data.merge$group) %>% addmargins()
levels(data.merge$orig.ident)
data.merge$orig.ident <- factor(data.merge$orig.ident,levels = c(paste0('Ctrl-',1:3),paste0('AD-',1:5)))
saveRDS(data.merge,file = 'data.merge.rds')
#step 2:QC
setwd('D:/科研工作站/单细胞工作站/H.B. Li/P-002/Integrated-3/Seurat')
ggplot(data.merge@meta.data,aes(x=orig.ident,y=nFeature_RNA,fill=orig.ident))+
  geom_violin()+
  theme_classic()+
  RotatedAxis()+theme(legend.position = 'none')
ggsave(filename = 'QC/vlnplot_nFeature_RNA.pdf')
ggplot(data.merge@meta.data,aes(x=orig.ident,y=nCount_RNA,fill=orig.ident))+
  geom_violin()+
  theme_classic()+
  RotatedAxis()+theme(legend.position = 'none')
ggsave(filename = 'QC/vlnplot_nCount_RNA.pdf')
ggplot(data.merge@meta.data,aes(x=orig.ident,y=percent.mt,fill=orig.ident))+
  geom_violin()+
  theme_classic()+
  RotatedAxis()+theme(legend.position = 'none')
ggsave(filename = 'QC/vlnplot_percent.mt.pdf')
data.merge <- subset(data.merge,subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 25)
table(data.merge$orig.ident,data.merge$group) %>% cbind() %>% addmargins()
dev.off()
VlnPlot(data.int,features = ,stack = T,group.by = 'orig.ident')+coord_flip()
VlnPlot(data.int,features = c('nCount_RNA','nFeature_RNA','percent.mt'),
        group.by = 'orig.ident',
        stack = T,
        flip = T)+
  theme(legend.position = 'none')+
  xlab('')+ylab('Value')
ggsave(filename = 'QC/merge.pdf',width = 6,height = 4)
#step 3:sketch a suset
data.merge <- NormalizeData(data.merge)
data.merge <- FindVariableFeatures(data.merge)
data.merge <- SketchData(
  object = data.merge,
  ncells = 5000,
  method = "LeverageScore",
  sketched.assay = "sketch")
DefaultAssay(data.merge) <- "sketch"
data.merge <- FindVariableFeatures(data.merge, verbose = F)
data.merge <- ScaleData(data.merge, verbose = F)
data.merge <- RunPCA(data.merge, verbose = F)
# integrate the datasets
library(future)
plan("multicore",works=12)
options(future.globals.maxSize = 4.8*1024^3 )
data.int <- IntegrateLayers(data.merge,
                            method = RPCAIntegration,
                            new.reduction = "integrated.rpca",
                            verbose = F)
# cluster the integrated data
ElbowPlot(data.int,ndims = 50)
ggsave('Clusters/ElbowPlot.pdf')
DimPlot(data.int,reduction = 'pca')
ggsave('Clusters/pca.pdf')
dev.off()
data.int <- FindNeighbors(data.int, reduction = 'integrated.rpca', dims = 1:30)
data.int <- FindClusters(data.int, resolution = seq(0.1,1,0.1))
data.int <- FindClusters(data.int, resolution = 0.4)
data.int <- RunUMAP(data.int, reduction = 'integrated.rpca', dims = 1:30,
                    return.model = T,
                    reduction.name = 'umap.rpca',
                    verbose = F)
colnames(data.int@meta.data)
DimPlot(data.int,group.by = 'sketch_snn_res.0.4',label = T)
ggsave(filename = 'Clusters/UMAP_0.4_clusters_sketch.pdf',width = 6,height = 5)
# Annotation
cellmarkers <- c('MZB1',
                 'CD79A',
                 'PECAM1',
                 'VWF',
                 'DCN',
                 'LUM',
                 'LYZ',
                 'HLA-DQA1',
                 'KLRD1',
                 'NCAM1',
                 'MKI67',
                 'TOP2A',
                 'ACTA2',
                 'TAGLN',
                 'CD3D',
                 'CD3E')
cellmarkers <- c('MZB1',
                 'PECAM1',
                 'DCN',
                 'LYZ',
                 'KLRD1',
                 'MKI67',
                 'ACTA2',
                 'CD3D')
colnames(md)
FeaturePlot(data.int,features = cellmarkers,reduction = 'umap.full',ncol = 4)
ggsave(filename = 'featureplot.pdf',width = 12,height = 5)
DotPlot(data.int,features = cellmarkers,group.by = 'sketch_snn_res.0.4')+coord_flip()
DotPlot(data.int,features = cellmarkers,group.by = 'celltype')+coord_flip()
data.int$celltype <- factor(data.int$celltype,levels = c('B',"EC","FB","Mφ",'NK','Proliferation',"SMC",'T','unclassified'))
data.int$celltype.full <- factor(data.int$celltype.full,levels = rev(c('B',"EC","FB","Mφ",'NK','Proliferation',"SMC",'T','unclassified')))
VlnPlot(data.int,
        features = cellmarkers,
        group.by = 'celltype',
        stack = T,
        flip = T,
        split.by = 'celltype',
        cols = brewer.pal(9,'Paired'))+
  theme(legend.position = 'none',
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())+
  xlab('')
ggsave(filename = 'Annotation/Vlnplot_0.4_clusters_sketch.pdf',height = 5,width = 8)
#FeaturePlot(data.int,features = cellmarkers,ncol = 3)
#ggsave(filename = 'Annotation/Featureplot_clusters_sketch.pdf',height = 12,width = 8)
Annotation <- readxl::read_xlsx('Annotation/Annotation.xlsx')
celltype <- Annotation$celltype %>% as.character()
subtype <- Annotation$subtype %>% as.character()
Idents(data.int) <- data.int$sketch_snn_res.0.4
levels(data.int)
names(celltype) <- levels(data.int)
data.int <- RenameIdents(data.int,celltype)
data.int <- AddMetaData(data.int,metadata = Idents(data.int),col.name = 'celltype')
Idents(data.int) <- data.int$sketch_snn_res.0.4
levels(data.int)
names(subtype) <- levels(data.int)
data.int <- RenameIdents(data.int,subtype)
data.int <- AddMetaData(data.int,metadata = Idents(data.int),col.name = 'subtype')
levels(data.int)
DimPlot(data.int, reduction = "umap.rpca",group='celltype',label = T)
ggsave(filename = 'Annotation/UMAP_0.4_celltype.sketch.pdf',width = 6,height = 5)
DimPlot(data.int, reduction = "umap.rpca",group='subtype')
ggsave(filename = 'Annotation/UMAP_0.4_subtype.nolable.sketch.pdf',width = 10,height = 5)
DimPlot(data.int, reduction = "umap.rpca",group='subtype',label = T)+NoLegend()
ggsave(filename = 'Annotation/UMAP_0.4_subtype_nolegend.sketch.pdf',width =5,height = 5)
#step 4:Project data to full cells
data.int <- ProjectIntegration(object = data.int, 
                               sketched.assay = "sketch", 
                               assay = "RNA",  
                               reduction = 'integrated.rpca')
data.int <- ProjectData(object = data.int, 
                        sketched.assay = "sketch", 
                        assay = "RNA",
                        sketched.reduction = 'integrated.rpca',
                        full.reduction = 'integrated.rpca.full', 
                        dims = 1:30, 
                        refdata = list(clusters.full = "sketch_snn_res.0.4",
                                       celltype.full='celltype',
                                       subtype.full='subtype'))
data.int <- RunUMAP(data.int, 
                    reduction = 'integrated.rpca', 
                    dims = 1:30, 
                    reduction.name = "umap.full",
                    reduction.key = "UMAP_full_")
DimPlot(data.int, reduction = "umap.full",group='clusters.full',label = T,label.size = 5,alpha = 0.8)+
  theme_test()+theme(plot.title = element_blank())+ylab('UMAP-2')+xlab('UMAP-1')+NoLegend()
ggsave(filename = 'Clusters/UMAP_0.4_clusters.full_nolegend.pdf',width = 5.5,height = 5)
DimPlot(data.int,group.by = 'clusters.full',label = T)+NoLegend()+xlab('UMAP-1')+ylab('UMAP-2')+theme(plot.title = element_blank())
DimPlot(data.int,group.by = 'celltype.full',label = T)+NoLegend()+xlab('UMAP-1')+ylab('UMAP-2')+theme(plot.title = element_blank())
ggsave(filename = 'dimplot_clusters.pdf',width = 5,height = 5)
saveRDS(data.int,file = 'data.int.rds')
#----------------------------------------Figure6: Hub genes expression patterns-------------------------------------------------------------------
#差异基因按细胞分组
deg_by_cell <- c()
for(i in 1:length(cell)){
  x <- FindMarkers(data.int,ident.1 = cell[i])
  x$change <- if_else(x$p_val_adj<0.05&abs(x$avg_log2FC)>1,if_else(x$avg_log2FC>1,'High','Low'),'Midium')
  x$change <- factor(x$change,levels = c('High','Low','Midium'))
  deg_by_cell[[i]] <- x}
names(data.int) <- cell
#hub gene表达
DefaultAssay(data.int) <- 'RNA'
Idents(data.int) <- data.int$celltype.full
hub <- c('HIF1A','HGF','ITGA5','HMOX1','ITGB3')
p1 <- DotPlot(data.int,features = hub)+theme_test()+xlab('')+ylab('')+theme()+coord_flip()
ggsave(filename = 'dotplot_expression_hub.pdf',width = 6.5,height = 5)
#差异基因按疾病分组
DefaultAssay(data.int) <- 'RNA'
data.int[['RNA']] <- JoinLayers(data.int[['RNA']])
data.int <- ScaleData(data.int)
cell <- levels(data.int)
deg_by_disease <- c()
for(i in 1:length(cell)){
  data <- subset(data.int,subset = celltype.full==cell[i])
  Idents(data) <- data$group
  x <- FindMarkers(data,ident.1 = 'AD',ident.2 = 'Ctrl')
  x$change <- if_else(x$p_val_adj<0.05&abs(x$avg_log2FC)>1,if_else(x$avg_log2FC>1,'Up','Down'),'Stable')
  x$change <- factor(x$change,levels = c('Up','Down','Stable'))
  deg_by_disease[[i]] <- x}
names(deg_by_disease) <- cell
saveRDS(deg_by_disease,file = 'deg_by_disease.rds')
#提取hub基因
degtf_by_disease <- c()
for(i in 1:length(deg_by_disease)){
  degtf_by_disease[[i]] <- deg_by_disease[[i]][tf,] %>% rownames_to_column('symbol')
}
names(degtf_by_disease) <- cell
writexl::write_xlsx(degtf_by_disease,path = 'deghub_by_disease.xlsx')
#可视化
degtf_by_disease <- readxl::read_xlsx('degtf_by_disease.xlsx')
colnames(degtf_by_disease)
degtf_by_disease$celltype <- factor(degtf_by_disease$celltype,levels = cell)
degtf_by_disease$symbol <- factor(degtf_by_disease$symbol,levels = tf)
p2 <- ggplot(degtf_by_disease,aes(x=celltype,y=symbol))+geom_point(aes(color=avg_log2FC,size=-log10(p_val_adj)))+
  scale_colour_gradient(low = 'lightgrey',high = 'red')+theme_test()+xlab('')+ylab('')
p1|p2
ggsave(filename = 'dotplot_expression_change.pdf',width = 14,height = 5)
#------------------------------------------Figure7 C: Analysis of HIF1A in macrophages -----------------------------------------------------------------------

#巨噬细胞差异基因，按疾病分组
deg.m <- deg_by_disease$Mφ
deg <- rownames(deg.m)[deg.m$change!=c('Stable')]
#可视化
deg.m1 <- deg.m %>% dplyr::filter(-log10(p_val_adj)<300)
ggplot(deg.m1,aes(x=avg_log2FC,y=-log10(p_val_adj),color=change))+
  geom_point()+
  scale_color_manual(values = c('#E16E1D','#3FA35C','gray'))+
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
  ylab('-log10(adj. P value)')+xlim(-5,5)
ggsave(filename = 'volcanoplot_deg.pdf',width = 6,height = 5)
#------------------------------------------Figure7 D: Analysis of HIF1A in macrophages -----------------------------------------------------------------------
#HIF1A的靶基因与巨噬细胞差异基因基因取交集
hif1a <- read_tsv(file = 'TRUST/HIF1A_targets.human.tsv')
hif1a <- hif1a$Target
intgene <- intersect(hif1a,deg)
deg_intgene <- deg.m[intgene,] %>% rownames_to_column('symbol')
writexl::write_xlsx(deg_intgene,path = 'deg_intgene_Mφ_HIF1A.xlsx')
#韦恩图
library(VennDiagram)
library(cowplot)
library(ggplotify)
vennplot <- venn.diagram(list('Mφ DEGs'=deg,
                              'HIF1A targeted genes'=hif1a),
                         filename = NULL,
                         cex=1.5,
                         col="transparent",
                         fill=c('#DD978F','skyblue'),
                         cat.cex=c(1,1),
                         cat.pos= c(0,0),
                         cat.dist=c(0.05,0.05),
                         scaled = F)
vennplot <- as.ggplot(as_grob(vennplot))
plot(vennplot)
ggsave(filename = 'VennDiagram-Mφ-HIF1A.pdf',width = 4,height = 4)

#--------------------------------------------------Figure8 D: HGF analysis-------------------------------------------------------------------------------
TKR <- c('NTRK1','MET','KDR','IGF1R','EGFR')
DotPlot(data.int,features = TKR)
ggsave(filename = 'dotplot_expression_HGFtar1')