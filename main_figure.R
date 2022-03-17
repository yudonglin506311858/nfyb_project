Sys.setenv(LANGUAGE = "en") #显示英文报错信息
gc()
memory.limit(9999999999)
set.seed(123)
rm(list = ls())  
options(stringsAsFactors = F)
library(Seurat)
library(tidyverse)
library(patchwork)
library(SCENIC)
library(plyr)
library(permute)
library(data.table)
library(SCopeLoomR)

setwd("D:/analysis/forpublication/human_bm_ubc/all")

YDL<-readRDS("D:/analysis/forpublication/human_bm_ubc/all/BM_all_SINGLET.RDS")


DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5)
DimPlot(YDL,reduction = "umap",label = TRUE,pt.size = 1.5)
DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5,group.by = "Phase")
DimPlot(YDL,reduction = "umap",label = TRUE,pt.size = 1.5,group.by = "Phase")
DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5,split.by = "orig.ident")
DimPlot(YDL,reduction = "umap",label = TRUE,pt.size = 1.5,split.by = "orig.ident")
DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5,group.by = "Phase",split.by = "orig.ident")
DimPlot(YDL,reduction = "umap",label = TRUE,pt.size = 1.5,group.by = "Phase",split.by = "orig.ident")




YDL.markers <- FindAllMarkers(YDL, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = 'roc')
# install.packages("magrittr") # package installations are only needed the first time you use it
# install.packages("dplyr")    # alternative installation of the %>%
library(magrittr) # needs to be run evYDL time you start R and want to use %>%
library(dplyr)    # alternatively, this also loads %>%
YDL.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
##存储marker

write.csv(YDL.markers,file="allmarker_human——cellytpe.csv")
#绘制分cluster的热图
top10 <- YDL.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_logFC)
#绘制marker在各个cluster的热图
pdf(file="tsneHeatmap_cellytpe.pdf",width=7,height=12)
DoHeatmap(object = YDL, features = top10$gene) + NoLegend()
DoHeatmap(subset(YDL, downsample = 100), features = top10$gene, size = 3)+ NoLegend()
dev.off()
table(YDL.markers$cluster)





library(clusterProfiler)
library("org.Hs.eg.db")
library(ggplot2)



#合并四个时期的BP，以热图展示


a <- read.csv("allmarker_human——cellytpe.csv")
b<-a[a$cluster=="ProE/BasoE","gene"]
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Hs.eg.db")
E1 <- eg[,2]
E1<-as.data.frame(E1)
colnames(E1)<-c("ProE/BasoE")
head(E1)
#E1<-c(E1)

a <- read.csv("allmarker_human——cellytpe.csv")
b<-a[a$cluster=="Early-PolyE","gene"]
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Hs.eg.db")
E2 <- eg[,2]
E2<-as.data.frame(E2)
colnames(E2)<-c("Early-PolyE")
head(E2)
#E2<-c(E2)

a <- read.csv("allmarker_human——cellytpe.csv")
b<-a[a$cluster=="Late-PolyE","gene"]
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Hs.eg.db")
E3 <- eg[,2]
E3<-as.data.frame(E3)
colnames(E3)<-c("Late-PolyE")
head(E3)
#E3<-c(E3)

a <- read.csv("allmarker_human——cellytpe.csv")
b<-a[a$cluster=="Early-OrthoE","gene"]
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Hs.eg.db")
E4 <- eg[,2]
E4<-as.data.frame(E4)
colnames(E4)<-c("Early-OrthoE")
head(E4)
#E4<-c(E4)

a <- read.csv("allmarker_human——cellytpe.csv")
b<-a[a$cluster=="Late-OrthoE","gene"]
eg = bitr(b, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Hs.eg.db")
E5 <- eg[,2]
E5<-as.data.frame(E5)
colnames(E5)<-c("Late-OrthoE")
head(E5)
#E4<-c(E4)


data<-list(`ProE/BasoE` =E1$`ProE/BasoE`, `Early-PolyE` = E2$`Early-PolyE`, `Late-PolyE`=E3$`Late-PolyE`,`Early-OrthoE` =E4$`Early-OrthoE`,`Late-OrthoE`=E5$`Late-OrthoE`)

lapply(data, head)


head(as.data.frame(ck))

ck <- compareCluster(geneCluster = data,OrgDb = org.Hs.eg.db, fun = "enrichGO", pvalueCutoff=0.05)
dotplot(ck, showCategory =20)
ego <- simplify(ck,cutoff=0.7,by="p.adjust",select_fun=min)
dotplot(ego, showCategory =15)
dotplot(ego, showCategory =20)+ RotatedAxis()

write.csv(ck,"ck_enrichGO.csv")

ekegg1<-setReadable(ck_enrichKEGG,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")


ck_enrichKEGG <- compareCluster(geneCluster = data, fun = "enrichKEGG",organism="hsa")
dotplot(ck_enrichKEGG, showCategory =25)
ck_enrichKEGG1<-setReadable(ck_enrichKEGG,OrgDb = org.Hs.eg.db,keyType = "SYMBOL")
write.csv(ck_enrichKEGG,"ck_enrichKEGG.csv")

#visualize the result using dotplot method.

pdf("多样本富集分析.pdf",width=10,height=8)
dotplot(ck, showCategory =10)

dotplot(ck, showCategory =50)


dotplot(ck, showCategory =20,split="ONTOLOGY") 
dev.off()








library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(magrittr)


DimPlot(YDL, reduction = "umap",pt.size = 1.5,label = T)
###Pseudotime monocle3
cds <- as.cell_data_set(YDL)
cds <- cluster_cells(cds)
head(pData(cds))
plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
plot_cells(cds, color_cells_by = "seurat_clusters", show_trajectory_graph = FALSE)
cds <- learn_graph(cds)
plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)


stem1<-rownames(pData(cds)[which(pData(cds)$ident %in% c('1')),])
cds <- order_cells(cds, root_cells = stem1)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, 
           label_branch_points = FALSE)


FeaturePlot(object = YDL, reduction = "umap",pt.size = 1.5,features = c("HBB","GATA1","GYPA","MALAT1"),cols = c("gray", "red"))#actin
FeaturePlot(object = YDL, reduction = "umap",pt.size = 1.5,features = c("HBB","GATA1","GYPA","KLF1"),cols = c("gray", "red"))#actin
FeaturePlot(object = YDL, reduction = "umap",pt.size = 1.5,features = c("XPO7"),cols = c("gray", "red"))#actin
FeaturePlot(object = YDL, reduction = "umap",pt.size = 1.5,features = c("BCL11A","NFYA","NFYB","NFYC"),cols = c("gray", "red"))#actin
FeaturePlot(object = YDL, reduction = "umap",pt.size = 1.5,features = c("XPO7"),cols = c("gray", "red"))#actin


VlnPlot(YDL,"HBB",slot = "data",pt.size =0)
VlnPlot(YDL,"GATA1",slot = "data",pt.size =0)
VlnPlot(YDL,"XPO7",slot = "data",pt.size =0)
VlnPlot(YDL,"NFYB",slot = "data",pt.size =0)
VlnPlot(YDL,"NFYC",slot = "data",pt.size =0)
VlnPlot(YDL,"TMCC2",slot = "data",pt.size =0)

#人工标记细胞类型
current.cluster.ids <- c(0, 1, 2, 3,4,5,6)

new.cluster.ids <- c(
  "Early-OrthoE",
  "ProE/BasoE",
  "Late-OrthoE",
  "Late-PolyE",
  "Late-OrthoE","Early-PolyE","Early-OrthoE")
names(new.cluster.ids) <- levels(YDL)
YDL <- RenameIdents(YDL, new.cluster.ids)
DimPlot(YDL, reduction = "umap", label = TRUE, pt.size = 1.5)

Idents(YDL)<-factor(Idents(YDL),levels =c("ProE/BasoE","Early-PolyE","Late-PolyE",
                                          "Early-OrthoE","Late-OrthoE"))

YDL$celltype<-Idents(YDL)




pdf("细胞类型指定-人.pdf")
DimPlot(YDL, reduction = "tsne", label = TRUE, pt.size = 1.5)
DimPlot(YDL,reduction = "tsne",label = TRUE,group.by="Phase",pt.size = 1.5)

DimPlot(YDL,reduction = "umap",label = TRUE,pt.size = 1.5)
DimPlot(YDL,reduction = "umap",label = TRUE,group.by="Phase",pt.size = 1.5)

dev.off()


YDL.AVERAGE<-AverageExpression(object = YDL,return.seurat=F)
YDL.AVERAGE<-as.data.frame(YDL.AVERAGE)

write.csv(YDL.AVERAGE,file = "AVERAGE_human.csv")
data<-YDL.AVERAGE
data<-as.matrix(data)
colnames(data)<-c("ProE/BasoE","Early-PolyE","Late-PolyE","Early-OrthoE","Late-OrthoE" )
head(data)
genes<-read.table("regulon.txt")
merge_tf<-c(genes$V1)
#挑选部分感兴趣
my.regulons <- merge_tf
#删掉所有列上都重复的
newdata<-data[c(my.regulons),]
newdata<-na.omit(newdata)
#colnames(newdata)<-c("ProE/BasoE","Early-PolyE","Late-PolyE","Early-OrthoE","Late-OrthoE" )
#低值为蓝色，高值为红色，中间值为白色：
#pdf("fig1g.pdf")
pheatmap(newdata,fontsize = 7,
         cluster_rows = F,border_color = NA,cluster_cols = F,color = colorRampPalette(colors = c("white","red","darkred"))(100))

pheatmap(newdata,scale = "row",fontsize = 7,
         cluster_rows = F,border_color = NA,cluster_cols = F,color = colorRampPalette(colors = c("navyblue","white","red"))(100))
pheatmap(newdata,scale = "row",fontsize = 7,
         cluster_rows = F,border_color = NA,cluster_cols = F,color = colorRampPalette(colors = c("blue","white","red"))(100))
pheatmap(newdata,scale = "row",fontsize = 7,
         cluster_rows = F,border_color = NA,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))





pheatmap(newdata,fontsize = 7,
         cluster_rows = T,border_color = NA,cluster_cols = T,color = colorRampPalette(colors = c("white","red","darkred"))(100))

pheatmap(newdata,scale = "row",fontsize = 7,
         cluster_rows = T,border_color = NA,cluster_cols = T,color = colorRampPalette(colors = c("navyblue","white","red"))(100))
pheatmap(newdata,scale = "row",fontsize = 7,
         cluster_rows = F,border_color = NA,cluster_cols = F,color = colorRampPalette(colors = c("blue","white","red"))(100))




# 加载包
library(Seurat)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
pheatmap(newdata,scale = "row",fontsize = 7,
         cluster_rows = F,border_color = NA,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))
pheatmap(newdata,scale = "row",fontsize = 7,
         cluster_rows = T,border_color = NA,cluster_cols = T,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))

pheatmap(newdata,scale = "row",fontsize = 7,clustering_method = "ward.D2",
         cluster_rows = T,border_color = NA,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))




pheatmap(newdata,scale = "row",fontsize = 6,filename = "new2.pdf",width = 10,height = 100,
         cluster_rows = F,border_color = NA,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))
dev.off()

pheatmap(data,fontsize = 7,
         cluster_rows = T,border_color = NA,cluster_cols = T,color = colorRampPalette(colors = c("white","red","darkred"))(100))

pheatmap(data,fontsize = 7,scale = "row",
         cluster_rows = T,border_color = NA,cluster_cols = T,color = colorRampPalette(colors = c("blue","white","red"))(100))

pheatmap(data,fontsize = 7,scale = "row",clustering_method = "ward.D",
         cluster_rows = T,border_color = NA,cluster_cols = T,color = colorRampPalette(colors = c("blue","white","red"))(100))
pheatmap(data,fontsize = 7,scale = "row",clustering_method = "ward.D",
         cluster_rows = T,border_color = NA,cluster_cols = F,color = colorRampPalette(colors = c("blue","white","red"))(100))
pheatmap(data,fontsize = 7,scale = "row",clustering_method = "ward.D",
         cluster_rows = T,border_color = NA,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))
pheatmap(data,fontsize = 7,scale = "row",clustering_method = "ward.D",
         cluster_rows = T,border_color = NA,cluster_cols = T,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))

pheatmap(data,fontsize = 7,scale = "row",clustering_method = "ward.D",show_rownames = T,
         cluster_rows = T,border_color = NA,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))

pheatmap(data,fontsize = 7,scale = "row",clustering_method = "ward.D2",show_rownames = T,
         cluster_rows = T,border_color = NA,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))









setwd("D:/analysis/forpublication/human_bm_ubc/all/SCENIC/int")

library(pheatmap)

library(RColorBrewer)
cellInfo<- read.table("Cell.Info.txt", sep = "\t", header = T, row.names = 1)

celltype = subset(cellInfo,select = 'CellType')
head(celltype)
#write.csv(celltype,"celltype.txt")
AUCmatrix<-read.table("AUCell.txt")
BINmatrix <- read.table("binary_mtx.txt", header = T)
AUCmatrix <- t(AUCmatrix)
BINmatrix <- t(BINmatrix)

colnames(AUCmatrix)<-celltype$CellType
library(gtools)
AUCmatrix<- AUCmatrix[,mixedorder(colnames(AUCmatrix))]
AUCmatrix[1:4,1:4]
dim(AUCmatrix)
table(colnames(AUCmatrix))

#计算每个cluster的score的平均值
cell_0<-apply(AUCmatrix[,colnames(AUCmatrix)%in%c("cell_0")],1,mean) #方差
cell_1<-apply(AUCmatrix[,colnames(AUCmatrix)%in%c("cell_1")],1,mean) #方差
cell_2<-apply(AUCmatrix[,colnames(AUCmatrix)%in%c("cell_2")],1,mean) #方差
cell_3<-apply(AUCmatrix[,colnames(AUCmatrix)%in%c("cell_3")],1,mean) #方差
cell_4<-apply(AUCmatrix[,colnames(AUCmatrix)%in%c("cell_4")],1,mean) #方差
cell_5<-apply(AUCmatrix[,colnames(AUCmatrix)%in%c("cell_5")],1,mean) #方差
cell_6<-apply(AUCmatrix[,colnames(AUCmatrix)%in%c("cell_6")],1,mean) #方差



data<-cbind(cell_0,cell_1,cell_2,cell_3,cell_4,cell_5,cell_6)

colnames(data)<-c(
  "Early-OrthoE",
  "ProE/BasoE",
  "Late-OrthoE",
  "Late-PolyE",
  "Late-OrthoE","Early-PolyE","Early-OrthoE")

head(data)
write.csv(data,"regulon_score.csv")



data<-data[which(rowSums(data) > 0),]#用R去除全是0的行


data<-read.csv("regulon_score.csv",row.names = 1)
data<-as.matrix(data)
colnames(data)<-c("ProE/BasoE","Early-PolyE","Late-PolyE","Early-OrthoE","Late-OrthoE" )
pheatmap(data,fontsize = 7,
         cluster_rows = T,border_color = NA,cluster_cols = T,color = colorRampPalette(colors = c("white","red","darkred"))(100))

pheatmap(data,fontsize = 7,scale = "row",
         cluster_rows = T,border_color = NA,cluster_cols = T,color = colorRampPalette(colors = c("blue","white","red"))(100))

pheatmap(data,fontsize = 7,scale = "row",clustering_method = "ward.D",
         cluster_rows = T,border_color = NA,cluster_cols = T,color = colorRampPalette(colors = c("blue","white","red"))(100))
pheatmap(data,fontsize = 7,scale = "row",clustering_method = "ward.D",
         cluster_rows = T,border_color = NA,cluster_cols = F,color = colorRampPalette(colors = c("blue","white","red"))(100))
pheatmap(data,fontsize = 7,scale = "row",clustering_method = "ward.D",
         cluster_rows = T,border_color = NA,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))
pheatmap(data,fontsize = 7,scale = "row",clustering_method = "ward.D",
         cluster_rows = T,border_color = NA,cluster_cols = T,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))

pheatmap(data,fontsize = 7,scale = "row",clustering_method = "ward.D",show_rownames = T,
         cluster_rows = T,border_color = NA,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))

pheatmap(data,fontsize = 7,scale = "row",clustering_method = "ward.D2",show_rownames = T,
         cluster_rows = T,border_color = NA,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))



