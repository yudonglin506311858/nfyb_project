sudo R
[sudo] password for ydl: 
R version 3.4.3

dir.create("results")
dir.create("results/paths")
#run this script in linux system and oracle java
source('/data/yudonglin/nopro/Cellrouter/cellrouter/CellRouter_Class.R')
libdir <- '/data/yudonglin/nopro/Cellrouter/cellrouter/CellRouter/'
set.seed(123)
library(dplyr)
library(plotrix)
matrix=read.table("projection.csv",sep=",",header=T,row.names=1)
colnames(matrix) <- c('UMAP1','UMAP2')
#rownames(matrix) = colnames (matrix)
rownames(matrix)=gsub("-1","",rownames(matrix))
ndata <- read.table('YDL.normalized_expression.txt',sep="\t",header=T,row.names=1)
genes <-as.vector(rownames(ndata))
map <- data.frame(id=rownames(ndata),symbol=genes,stringsAsFactors = FALSE)
ndata <- averageIds(ndata,map,'symbol')

#Remove genes with zero variance across all cells
var <- apply(ndata,1,var)
var <- var[which(var > 0)]
ndata <- ndata[names(var),]

### selecting genes to use as regulated along developmental trajectories.
#pca <- prcomp(t(ndata),scale=TRUE,center=TRUE)
#loadings <- pca$rotation
#num_pc <- 5
#quantile <- 0.975
#genes2use <- unique(as.vector(unlist(apply(loadings[,1:num_pc],2,function(x){names(x[which(abs(x) >= quantile(x,quantile))])}))))
genes2use=rownames(ndata)
ggrn <- buildGRN('Mm',ndata,genes2use,2,'results/GRN.R') #original 5
rownames(matrix) = colnames(ndata)



#下次可以直接load这个GRN文件
#ggrn <- get(load('results/GRN.R'))
### Subpopulation identification and gene signatures with CellRouter
cellrouter <- CellRouter(expdata=ndata,annotations=colnames(ndata))
cellrouter@rdimension <- matrix
pdf("kNN_network.pdf")
cellrouter <- findsubpopulations(cellrouter,90,'jaccard','results/kNN_network.gml')
dev.off()



df=read.table("YDL.meta.data.txt",header=T,row.names=1,sep="\t")
df$sample_id=rownames(df)
df=merge(cellrouter@sampTab,df,by="sample_id",all=T)
write.table(df,"YDL.meta.data.withSP.txt",sep="\t")

lengths(cellrouter@graph$subpopulation)
cellrouter <- diffexpr(cellrouter,column='population',pvalue = 0.05)
markers <- findmarkers(cellrouter)
write.table(markers,"results/YDL.markers.txt",sep="\t")
plotReducedDimension(cellrouter,5,5,filename='results/YDL.tSNE.pdf')
table(cellrouter@sampTab$population)
write.table(cellrouter@sampTab,"results/YDL.cellrouter_sampTab.txt",sep="\t")

######## Trajectory Detection using CellRouter ###
pdf("kNN_network_trajectory.pdf")
cellrouter <- createKNN(cellrouter,90,'jaccard','results/paths/kNN_network_trajectory.gml') #10 before this 90
dev.off()

filename <- "results/paths/cell_edge_weighted_network.txt"
write.table(cellrouter@graph$edges,file=filename,sep='\t',row.names=FALSE,col.names = FALSE,quote=FALSE) #input network
saveRDS(cellrouter,"cellrouter.RDS")

cellrouter<-readRDS("cellrouter.RDS")
##select starting subpopulation,all other subpopulations are targets
sources <- c('SP_9') #from SP_9 to SP_16
targets <- setdiff(as.vector(cellrouter@sampTab$population),sources)
methods <- c("euclidean","maximum","manhattan","canberra","binary",'graph') #graph for distances in KNN
cellrouter <- findpaths(cellrouter,libdir,paste(getwd(),'results/paths',sep='/'),method="graph")
ranks <- c('path_cost','path_flow','rank','length')
cellrouter <- processtrajectories(cellrouter,genes2use,path.rank=ranks[3],num.cells = 3,neighs = 1)
names <- unique(names(cellrouter@pathsinfo$distr))
clusters.show <- names
cellrouter <- correlationpseudotime(cellrouter,type='spearman')
cellrouter <- topgenes(cellrouter,0.85,0.15)
cellrouter <- smoothdynamics(cellrouter,names)
cellrouter <- clusterGenesPseudotime(cellrouter,10)
save(cellrouter,file='results/CellRouter_StemID_Processed.R')
saveRDS(cellrouter,"cellrouter_1.RDS")
cellrouter<- get(load('results/CellRouter_StemID_Processed.R'))

cellrouter<-readRDS("cellrouter_1.RDS")
##plot begins####
###positive and negative controls
p <- c('SP_9.SP_16') 
cellrouter@signatures$SP_1$subpopulation="SP_1"
cellrouter@signatures$SP_2$subpopulation="SP_2"
cellrouter@signatures$SP_3$subpopulation="SP_3"
cellrouter@signatures$SP_4$subpopulation="SP_4"
cellrouter@signatures$SP_5$subpopulation="SP_5"
cellrouter@signatures$SP_6$subpopulation="SP_6"
cellrouter@signatures$SP_7$subpopulation="SP_7"
cellrouter@signatures$SP_8$subpopulation="SP_8"
cellrouter@signatures$SP_9$subpopulation="SP_9"
cellrouter@signatures$SP_10$subpopulation="SP_10"
cellrouter@signatures$SP_11$subpopulation="SP_11"
cellrouter@signatures$SP_12$subpopulation="SP_12"
cellrouter@signatures$SP_13$subpopulation="SP_13"


data=rbind(cellrouter@signatures$SP_1,cellrouter@signatures$SP_2,cellrouter@signatures$SP_3,cellrouter@signatures$SP_4,cellrouter@signatures$SP_5,cellrouter@signatures$SP_6,cellrouter@signatures$SP_7,cellrouter@signatures$SP_8,cellrouter@signatures$SP_9,cellrouter@signatures$SP_10,cellrouter@signatures$SP_11,cellrouter@signatures$SP_12,cellrouter@signatures$SP_13)
data$gene_names=rownames(data)
write.table(data,"results/ydl.dif_genes.txt",sep="\t")
write.table(as.matrix(unlist(cellrouter@top.correlations$up)),"results/SP_6.SP_19_up_top_correlations.genes.txt",sep="\t")
write.table(as.matrix(unlist(cellrouter@top.correlations$down)),"results/SP_6.SP_19_down_top_correlations.genes.txt",sep="\t")

genlist=c("ALAS2","BNIP3L","SEC62","CA1","HISTIH4C","PTMA","TMCC2","ARL4A","BPGM","HIST1H4C","HMGB2","GPX1","KRT1","EIF1AY","MALAT1","MGST3","FTL","HBD")
plotPathHeatmap2(cellrouter,p,genelist,TRUE,2,2,10,10,paste('results/',"heatmap_along_trajectory__",sep=''))





## GRN score for selected transitions
tfs <- find_tfs(species = 'Hs')
save(tfs,file="results/tfs.R")

pdf("x.pdf")
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
dev.off()

pdf("results/grndynamics.cor.SP_6.SP_19.pdf")
grndynamics(cellrouter, tfs,p, 100)
dev.off()


scores <- x[[p]]$scores

pdf("results/grndynamics.cor.SP_6.SP_19.pdf")
plottrajectories(cellrouter, p, names(scores), rescale = TRUE, columns=1, width=20, height=10, filename='results/GRNscores_dynamics_curve.SP_6.SP_19.pdf') #trend lines
dev.off()


genelist <- c("Cebpd","Elk3","Sox4","Cd44")
pdf("tryall_2.pdf")
plotDRExpression(cellrouter,genelist,TRUE,2,10,10,paste('results/',p,"_some_genes_DRE.pdf",sep='')) #set less ploting genes,show marker exressing region 
plotheatmap(cellrouter, names(cellrouter@dynamics), 0.7, crows=5, ccols=1, width=10, height=10, filename="results/cor.path_heatmap.SP_10.SP_8.pdf") 
dotplot2(cellrouter,cellrouter@sampTab,markers_select,"population",1,logtransform=TRUE,15,15,filename="results/dotplot2.markers_select.SP_10.SP_8..pdf")
plotclusters(cellrouter, p,2,10,10, filename="results/plotclusters.SP_10.SP_8")
plotpaths(cellrouter, p, genelist, columns=7, width=23, height=10, file_prefix="results/plotpaths_1") #curve plus scatter means expression
plottrajectories(cellrouter,p,genelist,rescale = TRUE,columns=1,width=10,height=6,filename=paste('results/',p,'.genelist_1_dynamics_curve.SP_10.SP_8.pdf',sep=''))
dev.off()



transitions <- c('SP_6.SP_19','SP_20.SP_13','SP_14.SP_21','SP_6.SP_21')
pdf("grnscores.pdf")
x <- grnscores(cellrouter, tfs, transitions, direction='up', q.up=14, dir.targets='up', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_up')
dev.off()
p <- 'SP_9.SP_16'
scores <- x[[p]]$scores
m2 <- plottr(cellrouter, p, x[[p]]$scores, cluster=TRUE, 2, 2.5, 10, paste('results/', p, 'up_diff_dynamics.pdf',sep=''))





#transitions=names(cellrouter@pathsinfo$path)
#grntransition(cellrouter, tfs, transitions, dir.targets='up',q.up=0.95, q.down=0.05, columns=2, 50, 50, "results/grntransition")

cluster1_genes=names(cellrouter@clusters$SP_10.SP_8$clustering[which(cellrouter@clusters$SP_10.SP_8$clustering==1)])
cluster2_genes=names(cellrouter@clusters$SP_10.SP_8$clustering[which(cellrouter@clusters$SP_10.SP_8$clustering==2)])
cluster3_genes=names(cellrouter@clusters$SP_10.SP_8$clustering[which(cellrouter@clusters$SP_10.SP_8$clustering==3)])
cluster4_genes=names(cellrouter@clusters$SP_10.SP_8$clustering[which(cellrouter@clusters$SP_10.SP_8$clustering==4)])
cluster5_genes=names(cellrouter@clusters$SP_10.SP_8$clustering[which(cellrouter@clusters$SP_10.SP_8$clustering==5)])
cluster6_genes=names(cellrouter@clusters$SP_10.SP_8$clustering[which(cellrouter@clusters$SP_10.SP_8$clustering==6)])
cluster7_genes=names(cellrouter@clusters$SP_10.SP_8$clustering[which(cellrouter@clusters$SP_10.SP_8$clustering==7)])
cluster8_genes=names(cellrouter@clusters$SP_10.SP_8$clustering[which(cellrouter@clusters$SP_10.SP_8$clustering==8)])
cluster9_genes=names(cellrouter@clusters$SP_10.SP_8$clustering[which(cellrouter@clusters$SP_10.SP_8$clustering==9)])
cluster10_genes=names(cellrouter@clusters$SP_10.SP_8$clustering[which(cellrouter@clusters$SP_10.SP_8$clustering==10)])

genelist <- c("GYPA","EPOR","SNCA")
plotDRExpression(cellrouter,genelist,TRUE,2,10,10,paste('results/',p,"_some_genes_DRE.pdf",sep='')) #set less ploting genes,show marker exressing region 
plotheatmap(cellrouter, names(cellrouter@dynamics), 0.7, crows=5, ccols=1, width=10, height=10, filename="results/cor.path_heatmap.SP_10.SP_8.pdf") 
dotplot2(cellrouter,cellrouter@sampTab,markers_select,"population",1,logtransform=TRUE,15,15,filename="results/dotplot2.markers_select.SP_10.SP_8..pdf")
plotclusters(cellrouter, p,2,10,10, filename="results/plotclusters.SP_10.SP_8")
plotpaths(cellrouter, p, genelist, columns=7, width=23, height=10, file_prefix="results/plotpaths_1") #curve plus scatter means expression
plottrajectories(cellrouter,p,genelist,rescale = TRUE,columns=1,width=10,height=6,filename=paste('results/',p,'.genelist_1_dynamics_curve.SP_10.SP_8.pdf',sep=''))

#### Pathway enrichment analysis on selected trajectories
clustergenes=c(cluster1_genes,cluster2_genes,cluster3_genes,cluster4_genes,cluster5_genes,cluster6_genes,cluster7_genes,cluster8_genes,cluster9_genes,cluster10_genes)
write.table(clustergenes,"results/clustergenes.SP_10.SP_8.txt",sep="\t")
length(clustergenes)
ids <- read.table("results/ids.SP_10.SP_8.txt",header=T,sep="\t",stringsAsFactors=F) #cook right ids object of 10 cluster genes with geneIDannotation function
dim(ids)
colnames(ids)=c('entrezgene','external_gene_name','description','cytogenetic_location')
paths <- names(cellrouter@pathsinfo$path)
cellrouter <- pathwayenrichment(cellrouter,paths,cc=NULL,'human','org.Hs.eg.db',ids)
enr_UP <- pathwaycluster(cellrouter,cellrouter@pathwayenrichment$UP$GOBP,30,TRUE,30,20,'results/UP_GOBP.SP_10.SP_8.pdf')
enr_DOWN <- pathwaycluster(cellrouter,cellrouter@pathwayenrichment$DOWN$GOBP,30,TRUE,30,20,'results/DOWN_GOBP.SP_10.SP_8.pdf')
##regulatornetwork for several TFs
library('ggnetwork')
library('GGally')
library('geomnet')
library('network')
library('sna')
#regulators=c("FOXO3","RNF10","SNCA","ZNF737","PER1","IRF3")
regulators=c("TERF2IP","NFIX","RNF10","GTF2B","SNCA","PBX1")
regulatornetwork(x, regulators, 4, 2, 10, 15, 'results/regulator_networks.TERF2IP.pdf')







#run this script in linux system and oracle java
source('/data/yudonglin/nopro/Cellrouter/cellrouter/CellRouter_Class.R')
libdir <- '/data/yudonglin/nopro/Cellrouter/cellrouter/CellRouter/'
set.seed(123)
library(dplyr)
library(plotrix)
matrix=read.table("projection.csv",sep=",",header=T,row.names=1)
colnames(matrix) <- c('tSNE1','tSNE2')
#rownames(matrix) = colnames (matrix)
rownames(matrix)=gsub("-1","",rownames(matrix))
ndata <- read.table('YDL.normalized_expression.txt',sep="\t",header=T,row.names=1)
genes <-as.vector(rownames(ndata))
map <- data.frame(id=rownames(ndata),symbol=genes,stringsAsFactors = FALSE)
ndata <- averageIds(ndata,map,'symbol')

#Remove genes with zero variance across all cells
var <- apply(ndata,1,var)
var <- var[which(var > 0)]
ndata <- ndata[names(var),]

### selecting genes to use as regulated along developmental trajectories.
#pca <- prcomp(t(ndata),scale=TRUE,center=TRUE)
#loadings <- pca$rotation
#num_pc <- 5
#quantile <- 0.975
#genes2use <- unique(as.vector(unlist(apply(loadings[,1:num_pc],2,function(x){names(x[which(abs(x) >= quantile(x,quantile))])}))))
genes2use=rownames(ndata)
#ggrn <- buildGRN('Mm',ndata,genes2use,2,'results/GRN.R') #original 5
rownames(matrix) = colnames(ndata)

#saveRDS(cellrouter,"cellrouter.RDS")
#下次可以直接load这个GRN文件
ggrn <- get(load('results/GRN.R'))

cellrouter<-readRDS("cellrouter_1.RDS")

## GRN score for selected transitions
tfs <- find_tfs(species = 'Hs')
save(tfs,file="results/tfs.R")
tfs<-get(load('results/tfs.R'))
###positive and negative controls
p <- c('SP_9.SP_16') 
cellrouter@signatures$SP_1$subpopulation="SP_1"
cellrouter@signatures$SP_2$subpopulation="SP_2"
cellrouter@signatures$SP_3$subpopulation="SP_3"
cellrouter@signatures$SP_4$subpopulation="SP_4"
cellrouter@signatures$SP_5$subpopulation="SP_5"
cellrouter@signatures$SP_6$subpopulation="SP_6"
cellrouter@signatures$SP_7$subpopulation="SP_7"
cellrouter@signatures$SP_8$subpopulation="SP_8"
cellrouter@signatures$SP_9$subpopulation="SP_9"
cellrouter@signatures$SP_10$subpopulation="SP_10"
cellrouter@signatures$SP_11$subpopulation="SP_11"
cellrouter@signatures$SP_12$subpopulation="SP_12"
cellrouter@signatures$SP_13$subpopulation="SP_13"

p <- c('SP_9.SP_1')
pdf("SP_9.SP_1.pdf")
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
dev.off()

p <- c('SP_9.SP_2')
pdf("SP_9.SP_2.pdf")
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
dev.off()

p <- c('SP_9.SP_3')
pdf("SP_9.SP_3.pdf")
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
dev.off()

p <- c('SP_9.SP_4')
pdf("SP_9.SP_4.pdf")
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
dev.off()

p <- c('SP_9.SP_5')
pdf("SP_9.SP_5.pdf")
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
dev.off()

p <- c('SP_9.SP_6')
pdf("SP_9.SP_6.pdf")
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
dev.off()

p <- c('SP_9.SP_7')
pdf("SP_9.SP_7.pdf")
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
dev.off()

p <- c('SP_9.SP_8')
pdf("SP_9.SP_8.pdf")
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
dev.off()

p <- c('SP_9.SP_9')
pdf("SP_9.SP_9.pdf")
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
dev.off()

p <- c('SP_9.SP_10')
pdf("SP_9.SP_10.pdf")
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
dev.off()

p <- c('SP_9.SP_11')
pdf("SP_9.SP_11.pdf")
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
dev.off()

p <- c('SP_9.SP_12')
pdf("SP_9.SP_12.pdf")
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
dev.off()

p <- c('SP_9.SP_13')
pdf("SP_9.SP_13.pdf")
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap
dev.off()




p <- c('SP_9.SP_1')
pdf("grndynamics1.pdf")
grndynamics(cellrouter, tfs,p, 100)
dev.off()


transitions <- c('SP_9.SP_1','SP_9.SP_8','SP_9.SP_16')
pdf("grnscores_all.pdf")
x <- grnscores(cellrouter, tfs, transitions, direction='up', q.up=14, dir.targets='up', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_up')
dev.off()

transitions <- c('SP_9.SP_1','SP_9.SP_2')
pdf("grnscores_12.pdf")
x <- grnscores(cellrouter, tfs, transitions, direction='up', q.up=14, dir.targets='up', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_up')
x <- grnscores(cellrouter, tfs, transitions, direction='down', q.down=14, dir.targets='down', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_up')
dev.off()

transitions <- c('SP_9.SP_3')
pdf("grnscores_3.pdf")
x <- grnscores(cellrouter, tfs, transitions, direction='up', q.up=14, dir.targets='up', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_up')
dev.off()
transitions <- c('SP_9.SP_4','SP_9.SP_5','SP_9.SP_6')
pdf("grnscores_456.pdf")
x <- grnscores(cellrouter, tfs, transitions, direction='up', q.up=14, dir.targets='up', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_up')
x <- grnscores(cellrouter, tfs, transitions, direction='down', q.down=14, dir.targets='down', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_up')
dev.off()

transitions <- c('SP_9.SP_7','SP_9.SP_8')
pdf("grnscores_789.pdf")
x <- grnscores(cellrouter, tfs, transitions, direction='up', q.up=14, dir.targets='up', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_up')
x <- grnscores(cellrouter, tfs, transitions, direction='down', q.down=14, dir.targets='down', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_up')
dev.off()

transitions <- c('SP_9.SP_10','SP_9.SP_11','SP_9.SP_12','SP_9.SP_13')
pdf("grnscores_11121314.pdf")
x <- grnscores(cellrouter, tfs, transitions, direction='up', q.up=14, dir.targets='up', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_up')
x <- grnscores(cellrouter, tfs, transitions, direction='down', q.down=14, dir.targets='down', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_up')
dev.off()

transitions <- c('SP_9.SP_15','SP_9.SP_16','SP_9.SP_17')
pdf("grnscores_151617.pdf")
x <- grnscores(cellrouter, tfs, transitions, direction='up', q.up=14, dir.targets='up', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_up')
x <- grnscores(cellrouter, tfs, transitions, direction='down', q.down=14, dir.targets='down', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_up')
dev.off()



transitions <- c('SP_9.SP_1','SP_9.SP_2','SP_9.SP_3','SP_9.SP_4','SP_9.SP_5','SP_9.SP_6','SP_9.SP_7','SP_9.SP_8','SP_9.SP_9','SP_9.SP_10','SP_9.SP_11','SP_9.SP_12','SP_9.SP_13','SP_9.SP_14','SP_9.SP_15','SP_9.SP_16','SP_9.SP_17')
pdf("grnscores_sp.pdf")
x <- grnscores(cellrouter, tfs, transitions, direction='up', q.up=17, dir.targets='up', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_up')
x <- grnscores(cellrouter, tfs, transitions, direction='down', q.down=17, dir.targets='down', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_up')
dev.off()

transitions <- c('SP_9.SP_1','SP_9.SP_2','SP_9.SP_3','SP_9.SP_4')
pdf("grnscores_sp.pdf")
x <- grnscores(cellrouter, tfs, transitions, direction='up', q.up=14, dir.targets='up', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_up')
x <- grnscores(cellrouter, tfs, transitions, direction='down', q.down=14, dir.targets='down', columns=2, width=8, height=5, flip=FALSE, filename='results/lineage_regulators_score_up')
dev.off()

p <- c('SP_9.SP_1','SP_9.SP_2','SP_9.SP_4','SP_9.SP_5','SP_9.SP_6','SP_9.SP_7','SP_9.SP_8','SP_9.SP_10','SP_9.SP_11','SP_9.SP_12','SP_9.SP_13','SP_9.SP_14','SP_9.SP_15','SP_9.SP_16','SP_9.SP_17')
scores <- x[[p]]$scores
m2 <- plottr(cellrouter, p, x[[p]]$scores, cluster=TRUE, 2, 2.5, 10, paste('results/', p, 'up_diff_dynamics.pdf',sep=''))

p <- c('SP_9.SP_1')
scores <- x[[p]]$scores
pdf('SP_9_1.pdf')
m2 <- plottr(cellrouter, p, x[[p]]$scores, cluster=TRUE, 2, 2.5, 10, paste('results/', p, 'up_diff_dynamics.pdf',sep=''))
dev.off()
