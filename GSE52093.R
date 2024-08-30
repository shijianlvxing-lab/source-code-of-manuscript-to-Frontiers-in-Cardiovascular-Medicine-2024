#---------------------------------------Figure4 B: Boxplot of the hub genes in GSE52093----------------------------------------------------------------------
setwd('D:/科研工作站/生信工作站/My project/P003/GSE52093/deg')
exp <- readRDS("exp.rds")
md <- readRDS(file = 'md.rds')
deg <- readRDS(file = 'deg.rds')
#hub genes
#hub genes
hub <- c('HIF1A','HGF','ITGA5','HMOX1','ITGB3')
deg.hub <- deg[deg$symbol%in%hub,]
writexl::write_xlsx(deg.hub,path = 'deg.hub.xlsx')
#boxplot
exp.hub <- exp[hub,] %>% rownames_to_column('symbol')
exp.hub <- pivot_longer(exp.hub,cols = 2:13,names_to = 'sample',values_to = 'value')
exp.hub$group <- rep(md$Group,5)
ggplot(data = exp.hub,aes(x=symbol,y=value,fill=group))+geom_boxplot()+theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),legend.position = 'top')+xlab('')+ylab('Relative expression')+ggtitle('GSE52093')
ggsave(filename = 'boxplot-hubgene.pdf',width = 6,height = 3.5)
#---------------------------------------Figure4 C: ROC of the hub genes in GSE52093--------------------------------------------------------------------------
setwd('D:/科研工作站/生信工作站/My project/P003/GSE52093/roc')
library(pROC)
roc.data <- exp[hub,] %>% t() %>% as.data.frame()
group <- md$Group
roc <- list()
for(i in 1:ncol(roc.data)){
  roc[[i]] <- roc(group,roc.data[,i])
}
names(roc) <- colnames(roc.data)
roc.result <- data.frame(symbol=names(roc))
for(i in 1:1:ncol(roc.data)){
  roc.result[i,2] <- roc[[i]]$auc
}
colnames(roc.result)[2] <- 'auc'
roc.result <- roc.result %>% arrange(auc)
writexl::write_xlsx(roc.result,path = 'auc.xlsx')
saveRDS(roc,file = 'roc.rds')
#可视化，合并在一张图上
library(RColorBrewer)
col <- brewer.pal(7,'Set3')
pdf(file = 'roc.pdf',width = 6,height = 6)
plot(roc[[1]],col=col[1],lwd.ticks=3,grid=c(0.1,0.2),font=2,font.lab=2)
for(i in 2:length(roc)){lines(roc[[i]],col=col[i],lwd=3)}
lenged <- paste(roc.result$symbol,'(AUC=',round(roc.result$auc,2),')')
legend('bottomright',
       legend = lenged,
       col = col,
       lwd = 3)
dev.off()