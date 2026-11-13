library(Seurat)
library(reshape2)
library(dplyr)
library(ggthemes)
library(stringr)
library(SingleCellExperiment)
library(scater)
library(scran)
library(batchelor)

inrds = "inputFile.rds"
#输出路径
minCell = 10
kValue = 5

set.seed(2023)

pcadim = 1:10
computePCs = 50
nneighbors = 30
varGeneNum = 2000

resolution = 0.8

seuset= readRDS(inrds)
regressvar = seuset@commands$ScaleData.RNA$vars.to.regress

#将seuset按照样本拆分成list
seuset.list <- SplitObject(seuset, split.by = "orig.ident")

cellNum = lapply(seuset.list, ncol)
print(cellNum)
print(lapply(seuset.list,nrow))
#保留符合最小细胞数要求的seurat对象
seuset.list = seuset.list[cellNum>minCell]

#单个样本无法进行MNN分析，也没必要标准化
if (length(x = seuset.list) < 2) {
    stop("RDS file must contain multiple Seurat objects for MNN", call. = FALSE)
}
#合并样本
seuset = merge(
    x = seuset.list[[1]],
    y = seuset.list[2:length(x = seuset.list)]
)
seuset <- FindVariableFeatures(seuset, selection.method = "vst", nfeatures = varGeneNum)
scalegene = VariableFeatures(seuset)
seuset <- ScaleData(object = seuset, features = scalegene, vars.to.regress = regressvar)

objects.sce <- lapply(
    X = seuset.list,
    FUN = function(x, f) { 
      return(as.SingleCellExperiment(x = subset(x = x, features = VariableFeatures(seuset))))
    },
    f = VariableFeatures(seuset)
)

#开始进行MNN标化
print("RunMNN")
out <- do.call(fastMNN, c(objects.sce, list(k = kValue, d = computePCs)));

outcorrected = reducedDim(out)
rownames(outcorrected)=colnames(seuset)
colnames(outcorrected)=paste0("MNN_",1:ncol(outcorrected))
featureloading = as.matrix(rowData(out))
colnames(featureloading)=paste0("MNN_",1:ncol(featureloading))
#在seuset添加mnn结果
seuset[["mnn"]] <- CreateDimReducObject(
    embeddings = outcorrected,
    loadings = featureloading,
    assay = DefaultAssay(object = seuset),
    key = "mnn_"
)

print("RunUMAP")
seuset <- RunUMAP(seuset, dims = pcadim, n.neighbors = nneighbors, reduction.name = "umap3d", n.components = 3, reduction.key = "umap3d_",reduction="mnn")
seuset <- RunUMAP(seuset, dims = pcadim, n.neighbors = nneighbors,reduction="mnn")

print("RunTSNE")
seuset <- RunTSNE(object = seuset,dims = pcadim,reduction.name = "tsne3d",dim.embed=3, reduction.key = "tSNE3d_",check_duplicates = FALSE,reduction="mnn")
seuset <- RunTSNE(object = seuset,dims = pcadim,check_duplicates = FALSE,reduction="mnn")

methodValue = "GraphCluster"
seuset <- FindNeighbors(seuset, reduction = "mnn", dims = pcadim)
seuset <- FindClusters(seuset, resolution = resolution)

