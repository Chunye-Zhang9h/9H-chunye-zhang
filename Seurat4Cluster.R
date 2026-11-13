#Seurat4Cluster 请按需修改参数

options(stringAsFactors=FALSE)
#rds结果文件路径
rdsFile = "/inputfile.rds"
#输出结果前缀
prefix = "samplename"
#输出结果路径
outpath = "Clusteroutpath"
#选择一个聚类方法: graphbased,Kmeans,manual
clusterMethod = "graphbased"
#选择一个计算marker方法: MAST,wilcox,None(JustCluster)
MarkerGeneMethod = "wilcox"
#聚类方法选择kmean时的kmean值
Kmean = 10
#Logfc.Threshold卡值选择
threshold = 0.25
#最小pct卡值
minpct = 0.1
#聚类系数
resolution = 0.8
#输入指定聚类结果（Manual）,没有则为NULL
clusterFile = NULL
#算logfc的底数 2或e
base = "2"
OnlyPos = TRUE
returnMarker = TRUE

if(MarkerGeneMethod=="None"){returnMarker=FALSE}
if(!is.null(clusterFile)){ClusterMethod = "manual"}

if(base=="e"){
	base=exp(1)
}else if(base=="2"){
	base=2
}

library(SingleCellExperiment)
library(scater)
library(plyr)
library(reshape2)
library(Seurat)
library(mclust)
library(dplyr)
print("Start")
setwd(outpath)
seuset = readRDS(rdsFile)
assay = DefaultAssay(seuset)
if(ClusterMethod=="graphbased"){
	dir.create("GraphClust")
	if(!is.null(seuset@commands$RunUMAP)){
        dims = seuset@commands$RunUMAP@params$dims
    }else if(!is.null(seuset@commands$RunTSNE)){
		dims = seuset@commands$RunTSNE@params$dims
	}
	seuset <- FindNeighbors(seuset, reduction = "pca",dims = dims)
	seuset <- FindClusters(seuset, resolution = resolution)
	table = table(Idents(seuset),seuset@meta.data$orig.ident)
	print(table)
	table = dcast(data.frame(table),Var1~Var2)
	colnames(table)[1]="Cluster"
	write.table(table,file=paste0(prefix,"_GraphClust.Statistics.txt"),sep="\t",row.names=F,quote=F)
	data = data.frame(Cell = colnames(seuset),Cluster=Idents(seuset))
    write.table(data,file=paste0(prefix,"_GraphClust.Summary_Cell.txt"),sep="\t",row.names=F,quot=F)
	DefaultAssay(seuset) <- plotassay
	if(returnMarker){
        allmarkers <- FindAllMarkers(object = seuset, only.pos = onlypos, min.pct = minpct,logfc.threshold = threshold, test.use = MarkerGeneMethod,base=base)
		print(head(allmarkers))
		print(table(allmarkers$cluster))
		nrows=nrow(allmarkers)
	}else{
		nrows=0
	}
	if(nrows!=0){
		write.table(data.frame(allmarkers$gene,allmarkers[,1:6]),file=paste0(prefix,"_GraphClust.AllMarkerGenes.txt"),sep="\t",row.names=F,quote=F)
	}
	DefaultAssay(seuset) <- assay
	saveRDS(seuset, file = paste0("../",prefix,"_GraphClust.seuset.rds"))
}
if(ClusterMethod=="Kmeans"){
	i=Kmean
	if(!is.null(seuset@commands$RunUMAP)){
		dims = seuset@commands$RunUMAP@params$dims
	}else if(!is.null(seuset@commands$RunTSNE)){
        dims = seuset@commands$RunTSNE@params$dims
    }
    dir.create(paste0("Kmean",i))
    clusterid = as.character(kmeans(Embeddings(seuset,reduction="pca")[,dims],centers=i)$clust)
    order = sort(table(clusterid),decreasing = T)
    newclusters = 0:(length(table(clusterid))-1)
    names(newclusters) = names(order)
    clusterid = newclusters[clusterid]
	clusterid = as.factor(clusterid)
    Idents(seuset) <- clusterid
    
    table = table(Idents(seuset),seuset@meta.data$orig.ident)
    print(table)
    table = dcast(data.frame(table),Var1~Var2)
    colnames(table)[1]="Cluster"
    write.table(table,file=paste0(prefix,"_Kmean",i,".Statistics.txt"),sep="\t",row.names=F,quot=F)
    data = data.frame(Cell = colnames(seuset),Cluster=Idents(seuset))
    write.table(data,file=paste0(prefix,"_Kmean",i,".Summary_Cell.txt"),sep="\t",row.names=F,quot=F)
	DefaultAssay(seuset) <- plotassay
	if(returnMarker){
		allmarkers <- FindAllMarkers(object = seuset, only.pos = onlypos, min.pct = minpct,logfc.threshold = threshold, test.use = MarkerGeneMethod,base=base)
		print(head(allmarkers))
        print(table(allmarkers$cluster))
		nrows = nrow(allmarkers)
	}else{
		nrows = 0
	}
	if(nrows!=0){
		write.table(data.frame(allmarkers$gene,allmarkers[,1:6]),file=paste0(prefix,"_Kmean",i,".AllMarkerGenes.txt"),sep="\t",row.names=F,quote=F)

	}
	DefaultAssay(seuset) <- assay
	saveRDS(seuset, file = paste0("../",prefix,"_Kmean",i,".seuset.rds"))
}
if(ClusterMethod=="manual"){
	dir.create("Manual")
	manualcluster = read.table(args$clusterFile,header=T,sep="\t")
	clusterId = manualcluster[,2]
	names(clusterId) = manualcluster[,1]
	seuset = subset(seuset,cells=manualcluster[,1])
	seuset = SetIdent(object = seuset, cells = manualcluster[,1], clusterId)
	table = table(Idents(seuset),seuset@meta.data$orig.ident)
    print(table)
    table = dcast(data.frame(table),Var1~Var2)
    colnames(table)[1]="Cluster"
    write.table(table,file=paste0(prefix,"_Manual.Statistics.txt"),sep="\t",row.names=F,quot=F)
    data = data.frame(Cell = colnames(seuset),Cluster=Idents(seuset))
    write.table(data,file=paste0(prefix,"_Manual.Summary_Cell.txt"),sep="\t",row.names=F,quot=F)
	DefaultAssay(seuset) <- plotassay
	if(returnMarker){
        allmarkers <- FindAllMarkers(object = seuset, only.pos = onlypos, min.pct = minpct,logfc.threshold = threshold, test.use = MarkerGeneMethod,base=base)
        print(head(allmarkers))
        print(table(allmarkers$cluster))
		nrows = nrow(allmarkers)
	}else{
		nrows = 0
	}
	if(nrows!=0){
		write.table(data.frame(allmarkers$gene,allmarkers[,1:6]),file=paste0(prefix,"_Manual.AllMarkerGenes.txt"),sep="\t",row.names=F,quote=F)
	}
	DefaultAssay(seuset) <- assay
	saveRDS(seuset, file = paste0("../",prefix,"_Manual.seuset.rds"))
}
version = seuset@version
activeassay = seuset@active.assay
info = rbind(paste0("SeuratVersion:",packageVersion("Seurat")),
	paste0("SeuratObjectVersion:",packageVersion("SeuratObject")),
    paste0("DefaultAssay:",activeassay),
	paste0("PlotAssay:",plotassay),
	paste0("ClusterMethod:",ClusterMethod),
	paste0("MarkerGeneMethod:",MarkerGeneMethod)
)
if(ClusterMethod=="graphbased"){
	info = rbind(info,
	    paste0("pcadim:",paste(dims,collapse=",")),
		paste0("resolution:",resolution)
	)
}
if(ClusterMethod=="Kmeans"){
    info = rbind(info,
	    paste0("pcadim:",paste(dims,collapse=",")),
    	paste0("Kmean:",Kmean)
    )
}
if(MarkerGeneMethod!="None"){
	info = rbind(info,
        paste0("MinPct:",minpct),
		paste0("Logfc.Threshold:",threshold)
        )
}
write.table(info,file=paste0("../",prefix,".parameterInfo.txt"),sep="\t",quote=F,row.names=F,col.names=F)
