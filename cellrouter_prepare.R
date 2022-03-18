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


#先运行几个降维算法之后再读取文件
###prepare normalized expression for Cellrouter input
mydata<- FetchData(YDL,vars = c("UMAP_1","UMAP_2"))
write.csv(mydata,file = "projection.csv")
write.table(as.matrix(YDL@assays$RNA@counts),"YDL.normalized_expression.txt",sep="\t")
write.csv(colnames(YDL@assays$RNA@counts),"YDL.cell_names.csv")
write.csv(rownames(YDL@assays$RNA@counts),"YDL.gene_names.csv")
write.table(as.matrix(YDL@meta.data),"YDL.meta.data.txt")


