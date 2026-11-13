library(DoubletFinder)
library(Seurat)
library(tidyverse)

countsFile = "Sample.matrix"
Doubletrate = 0.1
pc.num = 1:10
pN = 0.25

allcounts = Read10X(data.dir=countsFile)
seuset <- CreateSeuratObject(counts = allcounts,min.cells = 3, min.features = 0)
seuset <- NormalizeData(seuset,normalization.method = "LogNormalize", scale.factor = 10000)
seuset <- FindVariableFeatures(seuset, selection.method = "vst", nfeatures = 2000)
scalegene = VariableFeatures(seuset)
seuset <- ScaleData(object = seuset,features = scalegene)
seuset <- RunPCA(seuset, verbose = F)
seuset <- RunUMAP(seuset, dims=pc.num)

seuset <- FindNeighbors(seuset, dims = pc.num) %>% FindClusters(resolution = resolution)
sweep.res.list <- paramSweep_v3(seuset, PCs = pc.num)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
homotypic.prop <- modelHomotypic(seuset$seurat_clusters)
nExp_poi <- round(Doubletrate*ncol(seuset))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

seuset <- doubletFinder_v3(seuset, PCs = pc.num, pN = pN, pK = pK_bcmvn, nExp = nExp_poi.adj, reuse.pANN = F)
