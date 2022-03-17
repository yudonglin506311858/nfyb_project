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

setwd("D:/analysis/forpublication/human_bm_ubc/all/SCENIC/int")


DimPlot(YDL,reduction = "umap",label = TRUE,pt.size = 1.5,split.by = "orig.ident",group.by="Phase")

DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5)
DimPlot(YDL,reduction = "umap",label = TRUE,pt.size = 1.5)
DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5,group.by = "Phase")
DimPlot(YDL,reduction = "umap",label = TRUE,pt.size = 1.5,group.by = "Phase")
DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5,split.by = "orig.ident")
DimPlot(YDL,reduction = "umap",label = TRUE,pt.size = 1.5,split.by = "orig.ident")
DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5,group.by = "Phase",split.by = "orig.ident")
DimPlot(YDL,reduction = "umap",label = TRUE,pt.size = 1.5,group.by = "Phase",split.by = "orig.ident")
#平均表达值
YDL.AVERAGE<-AverageExpression(object = YDL,return.seurat=F)
write.csv(YDL.AVERAGE,file = "AVERAGE_human.csv")



#读取细胞信息：
#人工标记细胞类型
current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6,7,8)
new.cluster.ids <- c(
  "cell_0",
  "cell_1",
  "cell_2",
  "cell_3",
  "cell_4",
  "cell_5",
  "cell_6")
names(new.cluster.ids) <- levels(YDL)
YDL <- RenameIdents(YDL, new.cluster.ids)
DimPlot(YDL, reduction = "tsne", label = TRUE, pt.size = 1.5)
DimPlot(YDL, reduction = "umap", label = TRUE, pt.size = 1.5)
x = Idents(object = YDL)
head(x)
y<-YDL@meta.data$seurat_clusters
head(y)
x<-as.data.frame(x)
head(x)
y<-as.data.frame(y)
head(y)
new<-cbind(x,y)
head(new)
colnames(new)<-c("CellType","ClusterID")
head(new)
write.table(new,"Cell.Info.txt",sep = "\t",quote = F)




cell.info <- read.table("Cell.Info.txt", sep = "\t", header = T, row.names = 1)
mydata_1<- FetchData(YDL,vars = c("tSNE_1","tSNE_2"))
mydata_2<- FetchData(YDL,vars = c("UMAP_1","UMAP_2"))
cell.info <- cbind(cell.info, mydata_1[rownames(cell.info), ])
cell.info <- cbind(cell.info, mydata_2[rownames(cell.info), ])
head(cell.info)
saveRDS(cell.info, "cell.info.rds")
cell.info <- readRDS("cell.info.rds")
## calculating the position of cluster labels
get_label_pos <- function(data, emb = "tSNE", group.by="ClusterID") {
  new.data <- data[, c(paste(emb, 1:2, sep = "_"), group.by)]
  colnames(new.data) <- c("x","y","cluster")
  clusters <- names(table(new.data$cluster))
  new.pos <- lapply(clusters, function(i) {
    tmp.data = subset(new.data, cluster == i)
    data.frame(
      x = median(tmp.data$x),
      y = median(tmp.data$y),
      label = i)
  })
  do.call(rbind, new.pos)
}
ggplot(cell.info, aes(tSNE_1, tSNE_2, color=as.character(ClusterID))) + 
  geom_point(size=.1) + 
  geom_text(inherit.aes = F, data = get_label_pos(cell.info, emb = "tSNE"), aes(x,y,label=label), size=3) + 
  theme_bw(base_size = 15) + 
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color="black")
  )
ggplot(cell.info, aes(UMAP_1, UMAP_2, color=as.character(ClusterID))) + 
  geom_point(size=.1) + 
  geom_text(inherit.aes = F, data = get_label_pos(cell.info, emb = "UMAP"), aes(x,y,label=label), size=3) + 
  theme_bw(base_size = 15) + 
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color="black")
  )
pdf(file="TSNE_umap_all.pdf",width=6.5,height=6)
TSNEPlot(object = YDL, pt.size = 2, label = TRUE)    #TSNE可视化
#另一个可视化的方法
DimPlot(YDL,reduction = "tsne",label = TRUE,pt.size = 1.5)
DimPlot(YDL,reduction = "umap",label = TRUE,pt.size = 1.5)
dev.off()

library(data.table)
library(pbapply)
library(plyr)
library(philentropy)
library(ggplot2)
library(ggrepel)
library(latex2exp)
setwd("./int")
#regulon在每个细胞中AUC值，最后得到一个以regulon为行细胞为列的矩阵。
regulonAUC<-readRDS("3.4_regulonAUC.Rds")
head(regulonAUC)#重点
AUCellThresholds<-readRDS("3.5_AUCellThresholds.Rds")
head(AUCellThresholds)
binaryRegulonActivity<-readRDS("4.1_binaryRegulonActivity.Rds")
head(binaryRegulonActivity[1:5,1:5])
binaryRegulonActivity_nonDupl<-readRDS("4.2_binaryRegulonActivity_nonDupl.Rds")
head(binaryRegulonActivity_nonDupl[1:5,1:5])
regulonSelections<-readRDS("4.3_regulonSelections.Rds")
head(regulonSelections)
write.table(regulonSelections$labels,"regulonSelections$labels.txt")
write.table(regulonSelections$all,"regulonSelections$all.txt")
write.table(regulonSelections$onePercent,"regulonSelections$onePercent.txt")
write.table(regulonSelections$corr,"regulonSelections$corr.txt")
write.table(regulonSelections$notCorr,"regulonSelections$notCorr.txt")
binaryRegulonOrder<-readRDS("4.4_binaryRegulonOrder.Rds")
head(binaryRegulonOrder)
write.table(binaryRegulonOrder,"binaryRegulonOrder.txt")
#regulonAUC<-as.data.frame(regulonAUC@assays$data)
regulonAUC<-as.data.frame(regulonAUC@assays@data$AUC)
#regulonAUC<-regulonAUC[,-1:-2]
regulonAUC<-t(regulonAUC)
head(regulonAUC[1:4,1:4])
dim(regulonAUC)
write.table(regulonAUC,"AUCell.txt")#更改(为_，删除)


getwd()
cell.info <- readRDS("cell.info.rds")
rasMat<-read.table("AUCell.txt")
head(rasMat[1:4,1:4])
cell.types <- names(table(cell.info$CellType))
ctMat <- lapply(cell.types, function(i) {
  as.numeric(cell.info$CellType == i)
})
ctMat <- do.call(cbind, ctMat)
colnames(ctMat) <- cell.types
rownames(ctMat) <- rownames(cell.info)
head(ctMat)
rssMat <- pblapply(colnames(rasMat), function(i) {
  sapply(colnames(ctMat), function(j) {
    1 - JSD(rbind(rasMat[, i], ctMat[, j]), unit = 'log2', est.prob = "empirical")
  })
})
rssMat <- do.call(rbind, rssMat)
rownames(rssMat) <- colnames(rasMat)
colnames(rssMat) <- colnames(ctMat)
saveRDS(rssMat, "rssMat.rds")
rssMat <- readRDS("rssMat.rds")

binaryRegulonActivity<-readRDS("4.1_binaryRegulonActivity.Rds")
head(binaryRegulonActivity[1:4,1:4])
binaryRegulonActivity<-t(binaryRegulonActivity)
head(binaryRegulonActivity[1:4,1:4])
write.table(binaryRegulonActivity,"binary_mtx.txt")

binMat <- read.table("binary_mtx.txt", header = T)
head(binMat[1:4,1:4])


PlotRegulonRank <- function(rssMat, cell.type, topn=5) {
  data <- data.frame(
    Regulons = 1:nrow(rssMat),
    RSS = sort(rssMat[, cell.type], decreasing = T),
    label = sub("(+)", "", names(sort(rssMat[, cell.type], decreasing = T)), fixed = T)
  )
  data$pt.col <- ifelse(data$Regulons <= topn, "#007D9B", "#BECEE3")
  data <- head(data, n=200)
  data.label <- head(data, n=topn)
  
  ggplot(data, aes(Regulons, RSS)) + 
    geom_point(size=3, color=data$pt.col) + 
    ggrepel::geom_text_repel(inherit.aes = FALSE, data = data.label, aes(Regulons, RSS, label=label), size=4) + 
    ggtitle(cell.type) + ylab("Specificity score") + 
    theme_bw(base_size = 12) + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color="black"),
          axis.ticks = element_line(color="black"),
          axis.text = element_text(color = "black"),
          plot.title = element_text(hjust = .5)
    )
}

data <- data.frame(
  Regulons = 1:nrow(rssMat),
  RSS = sort(rssMat[, "E0"], decreasing = T),
  label = sub("(+)", "", names(sort(rssMat[, "E0"], decreasing = T)), fixed = T)
)
write.csv(data,"E0_data.csv")

data <- data.frame(
  Regulons = 1:nrow(rssMat),
  RSS = sort(rssMat[, "E1"], decreasing = T),
  label = sub("(+)", "", names(sort(rssMat[, "E1"], decreasing = T)), fixed = T)
)
write.csv(data,"E1_data.csv")

data <- data.frame(
  Regulons = 1:nrow(rssMat),
  RSS = sort(rssMat[, "E2"], decreasing = T),
  label = sub("(+)", "", names(sort(rssMat[, "E2"], decreasing = T)), fixed = T)
)
write.csv(data,"E2_data.csv")


data <- data.frame(
  Regulons = 1:nrow(rssMat),
  RSS = sort(rssMat[, "E3"], decreasing = T),
  label = sub("(+)", "", names(sort(rssMat[, "E3"], decreasing = T)), fixed = T)
)
write.csv(data,"E3_data.csv")


data <- data.frame(
  Regulons = 1:nrow(rssMat),
  RSS = sort(rssMat[, "E3"], decreasing = T),
  label = sub("(+)", "", names(sort(rssMat[, "E3"], decreasing = T)), fixed = T)
)
write.csv(data,"E3_data.csv")

data <- data.frame(
  Regulons = 1:nrow(rssMat),
  RSS = sort(rssMat[, "E4"], decreasing = T),
  label = sub("(+)", "", names(sort(rssMat[, "E4"], decreasing = T)), fixed = T)
)
write.csv(data,"E4_data.csv")

data <- data.frame(
  Regulons = 1:nrow(rssMat),
  RSS = sort(rssMat[, "E5"], decreasing = T),
  label = sub("(+)", "", names(sort(rssMat[, "E5"], decreasing = T)), fixed = T)
)
write.csv(data,"E5_data.csv")



options(ggrepel.max.overlaps = Inf)
pdf("PlotRegulonRank_50.pdf")
PlotRegulonRank(rssMat, "E0",topn=10)
PlotRegulonRank(rssMat, "E1",topn=10)
PlotRegulonRank(rssMat, "E2",topn=10)
PlotRegulonRank(rssMat, "E3",topn=10)
PlotRegulonRank(rssMat, "E4",topn=10)
PlotRegulonRank(rssMat, "E5",topn=10)
dev.off()




cell.info <- readRDS("cell.info.rds")
cell.info <- cbind(cell.info, binMat[rownames(cell.info), ])
DimPlot <- function(cell.info, dim.1="UMAP_1", dim.2="UMAP_2", cell.type=NULL, regulon=NULL) {
  if (!is.null(cell.type)){
    data <- cell.info[, c(dim.1, dim.2, "CellType")]
    data$pt.col <- ifelse(data$CellType == cell.type, "red", "#DFDFDF")
    data$pt.size <- ifelse(data$CellType == cell.type, 0.2, 0.1)
    title <- paste0(cell.type)
    col.title <- "red"
  } else {
    data <- cell.info[, c(dim.1, dim.2, regulon)]
    data$pt.col <- ifelse(data[, regulon], "#006464", "#DFDFDF")
    data$pt.size <- ifelse(data[, regulon], 0.2, 0.1)
    title <- paste0("Regulon: ", regulon)
    col.title = "#006464"
  }
  ggplot(data, aes(get(dim.1), get(dim.2))) + 
    geom_point(size=data$pt.size, color=data$pt.col) + 
    theme_bw(base_size = 12) + 
    ggtitle("") + 
    xlab(dim.1) + ylab(dim.2) + 
    annotate("text",x=Inf,y=Inf,hjust=1.1,vjust=1.5,label=title,color=col.title,size=6) + 
    theme(legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color="black")
    )
}


pdf(file="regulon_tf_0.3_UMAP.pdf")
merge_tf<-read.csv("E0_data.csv")
merge_tf<-merge_tf[merge_tf$RSS>0.3,"label"]
for(i in merge_tf){ 
  p1<-DimPlot(cell.info, regulon = i)
  print(p1)
}
merge_tf<-read.csv("E1_data.csv")
merge_tf<-merge_tf[merge_tf$RSS>0.3,"label"]
for(i in merge_tf){ 
  p1<-DimPlot(cell.info, regulon = i)
  print(p1)
}
merge_tf<-read.csv("E2_data.csv")
merge_tf<-merge_tf[merge_tf$RSS>0.3,"label"]
for(i in merge_tf){ 
  p1<-DimPlot(cell.info, regulon = i)
  print(p1)
}
merge_tf<-read.csv("E3_data.csv")
merge_tf<-merge_tf[merge_tf$RSS>0.3,"label"]
for(i in merge_tf){ 
  p1<-DimPlot(cell.info, regulon = i)
  print(p1)
}
merge_tf<-read.csv("E4_data.csv")
merge_tf<-merge_tf[merge_tf$RSS>0.3,"label"]
for(i in merge_tf){ 
  p1<-DimPlot(cell.info, regulon = i)
  print(p1)
}
merge_tf<-read.csv("E5_data.csv")
merge_tf<-merge_tf[merge_tf$RSS>0.3,"label"]
for(i in merge_tf){ 
  p1<-DimPlot(cell.info, regulon = i)
  print(p1)
}
dev.off()



pdf("cluster.pdf")
DimPlot(cell.info, cell.type = "E0")
DimPlot(cell.info, cell.type = "E1")
DimPlot(cell.info, cell.type = "E2")
DimPlot(cell.info, cell.type = "E3")
DimPlot(cell.info, cell.type = "E4")
DimPlot(cell.info, cell.type = "E5")
dev.off()



write.csv(cell.info,"cell.info.csv")

DimPlot(cell.info, regulon = "FOXO3_67g")
dev.off()


figPlot <- function(cell.type, regulon) {
  p.list <- list(
    PlotRegulonRank(rssMat, cell.type),
    DimPlot(cell.info, cell.type = cell.type),
    DimPlot(cell.info, regulon = regulon)
  )
  cowplot::plot_grid(plotlist = p.list, ncol = 3, 
                     rel_widths = c(3,5,5))
}

pdf("regulon_merge.pdf",height=100,width=300)
figPlot(cell.type = "cell_0",regulon = "Phf8_extended_389g")
dev.off()

pdf(file="regulon_tf_results.pdf",height=10,width=30)
merge_tf<-read.csv("E0_data.csv")
merge_tf<-merge_tf[merge_tf$RSS>0.3,"label"]
for(i in merge_tf){ 
  p1<-figPlot(cell.type = "E0",regulon = i)
  print(p1)
}
merge_tf<-read.csv("E1_data.csv")
merge_tf<-merge_tf[merge_tf$RSS>0.3,"label"]
for(i in merge_tf){ 
  p1<-figPlot(cell.type = "E1",regulon = i)
  print(p1)
}
merge_tf<-read.csv("E2_data.csv")
merge_tf<-merge_tf[merge_tf$RSS>0.3,"label"]
for(i in merge_tf){ 
  p1<-figPlot(cell.type = "E2",regulon = i)
  print(p1)
}
merge_tf<-read.csv("E3_data.csv")
merge_tf<-merge_tf[merge_tf$RSS>0.3,"label"]
for(i in merge_tf){ 
  p1<-figPlot(cell.type = "E3",regulon = i)
  print(p1)
}
merge_tf<-read.csv("E4_data.csv")
merge_tf<-merge_tf[merge_tf$RSS>0.3,"label"]
for(i in merge_tf){ 
  p1<-figPlot(cell.type = "E4",regulon = i)
  print(p1)
}
merge_tf<-read.csv("E5_data.csv")
merge_tf<-merge_tf[merge_tf$RSS>0.3,"label"]
for(i in merge_tf){ 
  p1<-figPlot(cell.type = "E5",regulon = i)
  print(p1)
}
dev.off()

merge_tf<-read.csv("E0_data.csv")
label<-merge_tf$label
library(stringr)
pdf("REGULON_EXPRESSION.pdf")
for(i in label){ 
  feature<-str_sub(i,1,str_locate(i,"_")[1]-1)
  p1<-FeaturePlot(object = YDL, reduction = "tsne",pt.size = 1.5,features = feature,cols = c("gray", "red"))#actin
  print(p1)
}
dev.off()




merge_tf<-read.csv("cell_0_data.csv")
label<-merge_tf$label
# Seurat可视化SCENIC结果
# 
# 把SCENIC结果中最重要的regulonAUC矩阵导入Seurat，这样得到的可视化结果更容易与我们之前的分析联系起来。
scRNA<-YDL
##导入原始regulonAUC矩阵
AUCmatrix <- readRDS("3.4_regulonAUC.Rds")
AUCmatrix <- AUCmatrix@assays@data@listData$AUC
AUCmatrix <- data.frame(t(AUCmatrix), check.names=F)
RegulonName_AUC <- colnames(AUCmatrix)
RegulonName_AUC <- gsub(' \\(','_',RegulonName_AUC)
RegulonName_AUC <- gsub('\\)','',RegulonName_AUC)
colnames(AUCmatrix) <- RegulonName_AUC
scRNAauc <- AddMetaData(scRNA, AUCmatrix)
scRNAauc@assays$integrated <- NULL
#saveRDS(scRNAauc,'scRNAauc.rds')

##导入二进制regulonAUC矩阵
BINmatrix <- readRDS("4.1_binaryRegulonActivity.Rds")
BINmatrix <- data.frame(t(BINmatrix), check.names=F)
RegulonName_BIN <- colnames(BINmatrix)
RegulonName_BIN <- gsub(' \\(','_',RegulonName_BIN)
RegulonName_BIN <- gsub('\\)','',RegulonName_BIN)
colnames(BINmatrix) <- RegulonName_BIN
scRNAbin <- AddMetaData(scRNA, BINmatrix)
scRNAbin@assays$integrated <- NULL
#saveRDS(scRNAbin, 'scRNAbin.rds')

pdf("REGULONEXPRESSION_all.pdf", width=14 ,height=4)
for(i in label){ 
  try({p1 = FeaturePlot(scRNAbin, features=i, label=T, reduction = 'tsne',cols = c("gray", "red"))
  p2 = FeaturePlot(scRNAauc, features=i, label=T, reduction = 'tsne',cols = c("gray", "red"))
  feature<-str_sub(i,1,str_locate(i,"_")[1]-1)
  p3<-FeaturePlot(object = YDL, reduction = "tsne",pt.size = 1.5,features = feature,cols = c("gray", "red"))#actin
  plotc = p1|p2|p3
  print(plotc)})
}
dev.off()
pdf("REGULONEXPRESSION_all_umap2.pdf", width=14 ,height=4)

  for(i in label){ 
    try({
      p1 = FeaturePlot(scRNAbin, pt.size = 1.5,features=i, label=F, reduction = 'umap',cols = c("gray", "red"))
      p2 = FeaturePlot(scRNAauc,pt.size = 1.5, features=i, label=F, reduction = 'umap',cols = c("gray", "red"))
      feature<-str_sub(i,1,str_locate(i,"_")[1]-1)
      p3<-FeaturePlot(object = YDL, reduction = "umap",pt.size = 1.5,features = feature,cols = c("gray", "red"))#actin
      plotc = p1|p2|p3
      print(plotc)
    })
  }

dev.off()


rm(scRNAauc)
rm(scRNAbin)
rm(scRNA)

merge_tf<-read.csv("merge_tf.csv",header = F)#自己选择
label<-merge_tf$V1
pdf("REGULONEXPRESSION_tf_umap.pdf", width=14 ,height=4)
for(i in label){ 
  p1 = FeaturePlot(scRNAbin, features=i, label=T, reduction = 'umap',cols = c("gray", "red"))
  p2 = FeaturePlot(scRNAauc, features=i, label=T, reduction = 'umap',cols = c("gray", "red"))
  feature<-str_sub(i,1,str_locate(i,"_")[1]-1)
  p3<-FeaturePlot(object = YDL, reduction = "umap",pt.size = 1.5,features = feature,cols = c("gray", "red"))#actin
  plotc = p1|p2|p3
  print(plotc)
}
dev.off()


library(grid)
library(pbapply)
library(circlize)
library(ggsci)
library(plyr)
library(ggplot2)
library(ggrepel)
library(dendextend)
library(ComplexHeatmap)

dim(rasMat)
pccMat <- cor(rasMat)
CSI <- function(r1, r2) {
  delta <- pccMat[r1,r2]
  r.others <- setdiff(colnames(pccMat), c(r1,r2))
  N <- sum(pccMat[r1, r.others] < delta) + sum(pccMat[r2, r.others] < delta)
  M <- length(r.others) * 2
  return(N/M)
}
csiMat <- pblapply(rownames(pccMat), function(i) sapply(colnames(pccMat), function(j) CSI(i, j)))
csiMat <- do.call(rbind, csiMat)
rownames(csiMat) <- rownames(pccMat)
round(pccMat[1:10,1:10], 2)
round(csiMat[1:10,1:10], 2)
csiMat.binary <- matrix(as.numeric(csiMat >= 0.7), nrow = nrow(csiMat))
colnames(csiMat.binary) <- colnames(csiMat)
rownames(csiMat.binary) <- rownames(csiMat)
csiMat.binary[1:10,1:10]
saveRDS(csiMat, "csiMat.rds")
write.table(csiMat.binary, "csiMat.binary.txt", sep = "\t")
mat = readRDS("csiMat.rds")

h = 7
row_dend = as.dendrogram(hclust(dist(mat), method = "complete"))
clusters <- cutree(row_dend, h = h) # dendextend::cutree()
row_dend = color_branches(row_dend, h = h, col = pal_d3("category20")(20))
plot(row_dend)
col_range = c(0.7, 1)
col_fun <- colorRamp2(col_range, c("#FCF8DE", "#253177"))
ht <- Heatmap(
  matrix = mat,
  col = col_fun,
  name = "ht1",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_column_names = FALSE,
  show_row_names = FALSE,
  show_heatmap_legend = FALSE
)
lgd <- Legend(
  col_fun = col_fun, 
  title = "", 
  at = col_range, 
  labels = c("low", "high"), 
  direction = "horizontal",
  legend_width = unit(1, "in"),
  border = FALSE
)
draw(ht, heatmap_legend_list = list(lgd), heatmap_legend_side = c("bottom"))
decorate_heatmap_body("ht1", {
  tree = column_dend(ht)
  ind = cutree(as.hclust(tree), h = h)[order.dendrogram(tree)]
  first_index = function(l) which(l)[1]
  last_index = function(l) { x = which(l); x[length(x)] }
  clusters <- names(table(ind))
  x1 = sapply(clusters, function(x) first_index(ind == x)) - 1
  x2 = sapply(clusters, function(x) last_index(ind == x))
  x1 = x1/length(ind)
  x2 = x2/length(ind)
  grid.rect(x = x1, width = (x2 - x1), y = 1-x1, height = (x1 - x2), 
            hjust = 0, vjust = 0, default.units = "npc", 
            gp = gpar(fill=NA, col="#FCB800", lwd=3))
  grid.text(label = paste0("M",clusters),
            x = x2-length(clusters)/length(ind), y = 1-x1-(x2-x1)/2,
            default.units = "npc",
            hjust = 1, vjust = 0.5,
            gp = gpar(fontsize=12, fontface="bold"))
})
decorate_column_dend("ht1", {
  tree = column_dend(ht)
  ind = cutree(as.hclust(tree), h = h)[order.dendrogram(tree)]
  first_index = function(l) which(l)[1]
  last_index = function(l) { x = which(l); x[length(x)] }
  clusters <- names(table(ind))
  x1 = sapply(clusters, function(x) first_index(ind == x)) - 1
  x2 = sapply(clusters, function(x) last_index(ind == x))
  grid.rect(x = x1/length(ind), width = (x2 - x1)/length(ind), just = "left",
            default.units = "npc", gp = gpar(fill = pal_d3("category20")(20), alpha=.5, col = NA))
})

tree = column_dend(ht)
ind = cutree(as.hclust(tree), h = h)[order.dendrogram(tree)]
clusters <- names(table(ind))
regulon.clusters <- data.frame(regulon=names(ind), cluster=paste0("M",ind))
write.table(regulon.clusters, "regulon_clusters.txt", sep = "\t", quote = F, row.names = F)


regulon_M1<-regulon.clusters[regulon.clusters$cluster=="M1",]
write.table(regulon_M1, "regulon_M1_clusters.txt", sep = "\t", quote = F, row.names = F)
regulon_M2<-regulon.clusters[regulon.clusters$cluster=="M2",]
write.table(regulon_M2, "regulon_M2_clusters.txt", sep = "\t", quote = F, row.names = F)
regulon_M3<-regulon.clusters[regulon.clusters$cluster=="M3",]
write.table(regulon_M3, "regulon_M3_clusters.txt", sep = "\t", quote = F, row.names = F)
regulon_M4<-regulon.clusters[regulon.clusters$cluster=="M4",]
write.table(regulon_M4, "regulon_M4_clusters.txt", sep = "\t", quote = F, row.names = F)
regulon_M5<-regulon.clusters[regulon.clusters$cluster=="M5",]
write.table(regulon_M5, "regulon_M5_clusters.txt", sep = "\t", quote = F, row.names = F)
regulon_M6<-regulon.clusters[regulon.clusters$cluster=="M6",]
write.table(regulon_M6, "regulon_M6_clusters.txt", sep = "\t", quote = F, row.names = F)
regulon_M7<-regulon.clusters[regulon.clusters$cluster=="M7",]
write.table(regulon_M7, "regulon_M7_clusters.txt", sep = "\t", quote = F, row.names = F)
regulon_M8<-regulon.clusters[regulon.clusters$cluster=="M8",]
write.table(regulon_M8, "regulon_M8_clusters.txt", sep = "\t", quote = F, row.names = F)
regulon_M9<-regulon.clusters[regulon.clusters$cluster=="M9",]
write.table(regulon_M9, "regulon_M9_clusters.txt", sep = "\t", quote = F, row.names = F)


k = length(clusters)
cell.info <- readRDS("cell.info.rds")
moduleRasMat <- lapply(paste0("M",1:k), function(x){
  regulon.use <- subset(regulon.clusters, cluster == x)$regulon
  rowMeans(rasMat[, regulon.use])
})
names(moduleRasMat) <- paste0("M",1:k)
moduleRasMat <- do.call(cbind, moduleRasMat)
cell.info <- cbind(cell.info, moduleRasMat[rownames(cell.info), ])
p.list <- lapply(paste0("M",1:k), function(module){
  data.use <- cell.info
  expression.color <- c("darkblue", "lightblue", "green", "yellow", "red")
  max.val <- quantile(data.use[, module], 0.99)
  low.val <- quantile(data.use[, module], 0.1)
  data.use[, module] <- ifelse(data.use[, module] > max.val, max.val, data.use[, module])
  ggplot(data.use, aes(tSNE_1, tSNE_2, color=get(module))) + 
    geom_point(size=0.05) + 
    theme_bw(base_size = 15) + 
    ggtitle(module) + 
    scale_color_gradientn(name = NULL, colors = expression.color) + 
    theme(legend.position = "right",
          legend.title = element_blank(),
          plot.title = element_text(hjust = .5, face = "bold", size = 20),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color="black")
    )
})
pdf("module_RESULT.pdf",width = 20,height = 20)
cowplot::plot_grid(plotlist = p.list, ncol = 3)
dev.off()
CellType.info <- read.table("Cell.Info.txt", sep = "\t", header = T)
## group module score by clusterID
clusters <- names(table(cell.info$ClusterID))
moduleScoreInCellType <- lapply(clusters, function(x) {
  cells.use <- rownames(subset(cell.info, ClusterID == x))
  colMeans(cell.info[cells.use, paste0("M",1:k)])
})
names(moduleScoreInCellType) <- clusters
moduleScoreInCellType <- do.call(rbind, moduleScoreInCellType)
moduleScoreInCellType <- as.data.frame(moduleScoreInCellType)
moduleScoreInCellType$CellType <- mapvalues(
  x = rownames(moduleScoreInCellType),
  from = CellType.info$Cluster,
  to = CellType.info$CellType
)
plotCellTypeRank <- function(data, module, topn=5){
  data.use <- data
  data.use <- data.use[order(data.use[, module], decreasing = TRUE), ]
  data.use$Rank <- 1:nrow(data.use)
  data.use$pt.col <- ifelse(data.use$Rank <= topn, "#007D9B", "#BECEE3")
  data.label <- head(data.use, n = topn)
  data.label$delta <- c(Inf, abs(diff(data.label[, module])))
  ggplot(data.use, aes(Rank, get(module))) + 
    geom_point(size=3, color=data.use$pt.col) + 
    geom_text_repel(inherit.aes = FALSE, data = data.label, aes(Rank, get(module), label=CellType), size=4, max.iter = 2e4) + 
    ggtitle(module) + 
    ylab("Regulon activity score") + 
    xlab("Cell type") + 
    theme_bw(base_size = 12) + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color="black"),
          axis.ticks = element_line(color="black"),
          axis.text = element_text(color = "black"),
          plot.title = element_text(hjust = .5, face = "bold")
    )
}
p.list <- lapply(paste0("M",1:k), function(x) plotCellTypeRank(moduleScoreInCellType, module = x, topn = 10))
cowplot::plot_grid(plotlist = p.list, ncol = 4)

pdf("module.pdf",height=100,width=20)
cowplot::plot_grid(plotlist = p.list, ncol = 4)
dev.off()


pdf("new.pdf",height=100,width =20)
cowplot::plot_grid(plotlist = p.list, ncol = 1)
dev.off()
savehistory("regulon_code.txt")




#pheatmap可视化SCENIC结果

library(pheatmap)
cellInfo<- read.table("Cell.Info.txt", sep = "\t", header = T, row.names = 1)

celltype = subset(cellInfo,select = 'CellType')
AUCmatrix<-read.table("AUCell.txt")
BINmatrix <- read.table("binary_mtx.txt", header = T)
AUCmatrix <- t(AUCmatrix)
BINmatrix <- t(BINmatrix)

merge_tf<-read.csv("E0_data.csv",header = F)#自己选择

merge_tf<-c(merge_tf$V1)
#挑选部分感兴趣的regulons
my.regulons <- merge_tf
myAUCmatrix <- AUCmatrix[rownames(AUCmatrix)%in%my.regulons,]
myBINmatrix <- BINmatrix[rownames(BINmatrix)%in%my.regulons,]
#使用regulon原始AUC值绘制热图
pheatmap(myAUCmatrix,fontsize = 5, show_colnames=F, annotation_col=celltype,cluster_rows = T,cluster_cols = T)#, show_colnames=F,
myAUCmatrix_1<-scale(myAUCmatrix,scale = T)
pheatmap(myAUCmatrix,scale = "row",fontsize = 5, show_colnames=F, annotation_col=celltype,cluster_rows = T,cluster_cols = T,color = colorRampPalette(colors = c("blue","red","red"))(100))#, show_colnames=F,
#filename = 'myAUCmatrix_heatmap.png',
#width = 6, height = 5)
#使用regulon二进制AUC值绘制热图
pheatmap(myBINmatrix, show_colnames=F,fontsize = 5, annotation_col=celltype,cluster_rows = T,
         color = colorRampPalette(colors = c("white","black"))(100))#,
filename = 'myBINmatrix_heatmap.pdf',
width = 6, height = 5)
dev.off()




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
library(Matrix)
library(cowplot)
library(harmony)
set.seed(123)





pdf("gene.pdf")
merge_tf<-read.csv("regulon_exp.csv",header = F)
merge_tf<-read.csv("regulon_gene.csv",header = F)
label<-merge_tf$V1
FeaturePlot(object = YDL, reduction = "umap",pt.size = 1.5,features = c(label),cols = c("gray", "red"))#actin
merge_tf<-read.csv("regulon_exp.csv",header = F)
label<-merge_tf$V1
FeaturePlot(YDLbin, features=c(label), label=T, reduction = 'tsne',cols = c("gray", "red"))
FeaturePlot(YDLauc, features=c(label), label=T, reduction = 'tsne',cols = c("gray", "red"))
dev.off()
# Seurat可视化SCENIC结果
# 
# 把SCENIC结果中最重要的regulonAUC矩阵导入Seurat，这样得到的可视化结果更容易与我们之前的分析联系起来。

##导入原始regulonAUC矩阵
AUCmatrix <- readRDS("3.4_regulonAUC.Rds")
AUCmatrix <- AUCmatrix@assays@data@listData$AUC
AUCmatrix <- data.frame(t(AUCmatrix), check.names=F)
RegulonName_AUC <- colnames(AUCmatrix)
RegulonName_AUC <- gsub(' \\(','_',RegulonName_AUC)
RegulonName_AUC <- gsub('\\)','',RegulonName_AUC)
colnames(AUCmatrix) <- RegulonName_AUC
YDLauc <- AddMetaData(YDL, AUCmatrix)
YDLauc@assays$integrated <- NULL
#saveRDS(YDLauc,'YDLauc.rds')

##导入二进制regulonAUC矩阵
BINmatrix <- readRDS("4.1_binaryRegulonActivity.Rds")
BINmatrix <- data.frame(t(BINmatrix), check.names=F)
RegulonName_BIN <- colnames(BINmatrix)
RegulonName_BIN <- gsub(' \\(','_',RegulonName_BIN)
RegulonName_BIN <- gsub('\\)','',RegulonName_BIN)
colnames(BINmatrix) <- RegulonName_BIN
YDLbin <- AddMetaData(YDL, BINmatrix)
YDLbin@assays$integrated <- NULL
#saveRDS(YDLbin, 'YDLbin.rds')

pdf("REGULONEXPRESSION_all.pdf", width=14 ,height=4)
for(i in label){ 
  p1 = FeaturePlot(YDLbin, features=i, label=T, reduction = 'tsne',cols = c("gray", "red"))
  p2 = FeaturePlot(YDLauc, features=i, label=T, reduction = 'tsne',cols = c("gray", "red"))
  feature<-str_sub(i,1,str_locate(i,"_")[1]-1)
  p3<-FeaturePlot(object = YDL, reduction = "tsne",pt.size = 1.5,features = feature,cols = c("gray", "red"))#actin
  plotc = p1|p2|p3
  print(plotc)
}
dev.off()
#pheatmap可视化SCENIC结果

library(pheatmap)
cellInfo<- read.table("Cell.Info.txt", sep = "\t", header = T, row.names = 1)

celltype = subset(cellInfo,select = 'CellType')
head(celltype)
write.csv(celltype,"celltype.txt")
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
E0<-apply(AUCmatrix[,colnames(AUCmatrix)%in%c("E0")],1,mean) #方差
E1<-apply(AUCmatrix[,colnames(AUCmatrix)%in%c("E1")],1,mean) #方差
E2<-apply(AUCmatrix[,colnames(AUCmatrix)%in%c("E2")],1,mean) #方差
E3<-apply(AUCmatrix[,colnames(AUCmatrix)%in%c("E3")],1,mean) #方差
E4<-apply(AUCmatrix[,colnames(AUCmatrix)%in%c("E4")],1,mean) #方差
E5<-apply(AUCmatrix[,colnames(AUCmatrix)%in%c("E5")],1,mean) #方差




colnames(BINmatrix)<-celltype$CellType
library(gtools)
BINmatrix<- BINmatrix[,mixedorder(colnames(BINmatrix))]
BINmatrix[1:4,1:4]
dim(BINmatrix)
table(colnames(BINmatrix))
#计算每个cluster的score的平均值
E0<-apply(BINmatrix[,colnames(BINmatrix)%in%c("E0")],1,mean) #方差
E1<-apply(BINmatrix[,colnames(BINmatrix)%in%c("E1")],1,mean) #方差
E2<-apply(BINmatrix[,colnames(BINmatrix)%in%c("E2")],1,mean) #方差
E3<-apply(BINmatrix[,colnames(BINmatrix)%in%c("E3")],1,mean) #方差
E4<-apply(BINmatrix[,colnames(BINmatrix)%in%c("E4")],1,mean) #方差
E5<-apply(BINmatrix[,colnames(BINmatrix)%in%c("E5")],1,mean) #方差


data<-cbind(E0,E1,E2,E3,E4,E5)
#只选择晚期的cluster:0,2,5,7,8
data<-data[,c(1,8,3,6,9)]

#展示这些regulon的表达水平
regulon_expression<-read.table("regulon_expression.txt")
features<-regulon_expression$V1
#绘制marker在各个cluster的热图
pdf(file="tsneHeatmap_all_regulon.pdf",width=12,height=9)
DoHeatmap(object = YDL, features = features) + NoLegend()
DoHeatmap(subset(YDL, downsample = 100), features = features, size = 3)+ NoLegend()
dev.off()

#merge_tf<-read.csv("merge_tf.csv",header = F)#自己选择
#merge_tf<-read.csv("merge_tf1.csv",header = F)#自己选择
#merge_tf<-read.csv("merge_tf2.csv",header = F)#自己选择
merge_tf<-read.csv("E0_data.csv",header = T)#自己选择

merge_tf<-c(merge_tf$label)
#挑选部分感兴趣的regulons
my.regulons <- merge_tf
#删掉所有列上都重复的
newdata<-data[c(my.regulons),]
#低值为蓝色，高值为红色，中间值为白色：
#pdf("fig1g.pdf")
pheatmap(newdata,fontsize = 7,
         cluster_rows = F,border_color = NA,cluster_cols = F,color = colorRampPalette(colors = c("white","red","darkred"))(100))

pheatmap(newdata,scale = "row",fontsize = 7,
         cluster_rows = F,border_color = NA,cluster_cols = F,color = colorRampPalette(colors = c("navyblue","white","red"))(100))
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
pheatmap(newdata,scale = "row",fontsize = 7,show_rownames = F,
         cluster_rows = T,border_color = NA,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))



pheatmap(newdata,scale = "row",fontsize = 6,filename = "new2.pdf",
         cluster_rows = F,border_color = NA,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))
dev.off()


#pheatmap可视化SCENIC结果

library(pheatmap)
cellInfo<- read.table("Cell.Info.txt", sep = "\t", header = T, row.names = 1)

celltype = subset(cellInfo,select = 'CellType')
head(celltype)
write.csv(celltype,"celltype.txt")
AUCmatrix<-read.table("AUCell.txt")
BINmatrix <- read.table("binary_mtx.txt", header = T)
AUCmatrix <- t(AUCmatrix)
BINmatrix <- t(BINmatrix)
merge_tf<-read.csv("merge_tf.csv",header = F)#自己选择
merge_tf<-c(merge_tf$V1)
#挑选部分感兴趣的regulons
my.regulons <- merge_tf
myAUCmatrix <- AUCmatrix[rownames(AUCmatrix)%in%my.regulons,]
myBINmatrix <- BINmatrix[rownames(BINmatrix)%in%my.regulons,]
#使用regulon原始AUC值绘制热图
pheatmap(myAUCmatrix,fontsize = 5, show_colnames=F, annotation_col=celltype,cluster_rows = T,cluster_cols = T)#, show_colnames=F,
myAUCmatrix_1<-scale(myAUCmatrix,scale = T)
pheatmap(myAUCmatrix,scale = "row",fontsize = 5, show_colnames=F, annotation_col=celltype,cluster_rows = T,cluster_cols = T,color = colorRampPalette(colors = c("blue","red","red"))(100))#, show_colnames=F,
#filename = 'myAUCmatrix_heatmap.png',
#width = 6, height = 5)
#使用regulon二进制AUC值绘制热图
pheatmap(myBINmatrix, show_colnames=F,fontsize = 5, annotation_col=celltype,cluster_rows = T,
         color = colorRampPalette(colors = c("white","black"))(100),
         filename = 'myBINmatrix_heatmap.pdf',
         width = 6, height = 5)
dev.off()
