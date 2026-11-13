library(monocle)
library(Seurat)
library(reshape2)
library(dplyr)

rdsFile = "./seuset.rds";
outpath = "./";
geneFile = "./gene.txt";
celllist = "./celllist.txt";

setwd(outpath)
seuset = readRDS(rdsFile)

seuset <- AddMetaData(object = seuset, metadata = seuset@active.ident, col.name = "Cluster")
seuCDSAll = importCDS(seuset, import_all = TRUE)
if(!is.null(geneFile)){
	genelist = unique(read.table(geneFile,sep="\t",header=T)[,1])
}else{
	genelist = VariableFeatures(seuset)
}

seuCDSAll <- setOrderingFilter(seuCDSAll, genelist)
seuCDSAll <- estimateSizeFactors(seuCDSAll)
seuCDS <- reduceDimension(seuCDSAll, max_components = 2, method = "DDRTree")
seuCDS <- orderCells(seuCDS)

plot_cell_trajectory(seuCDS, color_by = "State")+theme(axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
plot_cell_trajectory(seuCDS, color_by = "Cluster")+theme(axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
plot_cell_trajectory(seuCDS, color_by = "Pseudotime")+theme(axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))

