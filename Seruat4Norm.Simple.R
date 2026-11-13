#利用Seurat对cellranger结果进行分析标化

countsFile = "MatrixDir"
#选择matrix结果中 gene表读第一列还是第二列
genecolumn = 2
#输出名称前缀
prefix = "SampleName"
#输出路径
outpath = "/Outpath"
#是否添加样本信息Samplefile
samplelist = "sample.txt"
#默认200每个细胞最少有200个基因表达，否则过滤掉这个细胞。
cellMinGene = 200
#默认10000表示一个细胞表达基因超过5000有可能是双细胞，被过滤掉
cellMaxGene = 10000
#默认3，表示一个基因至少要在三个细胞里有表达量，否则过滤掉这个基因。
geneExpMinCell = 3
#默认0.3表示线粒体基因超过20%有可能是一个死细胞，会被过滤掉
MaxMTPercent = 0.3

#选择要删除的gene:MTGene;
removeGene = "MTGene"
#数据中心化处理的算法 Seurat提供的3种不同算法 默认linear 运算速度最快 poisson negbinom
ScaleModel = "linear"
#标化方法： scale或者 sctransform(SCTransform)
NormType = "scale"         

regressOut = "nCount_RNA,percent.mito"
#选择要进行的额外分析：多个用逗号隔开

runFast=TRUE
if(NormType=="scale"){
	pcadim = 1:10
	computePCs = 50
	nneighbors = 30
	varGeneNum = 2000		
}else if(NormType=="sctransform"){
	pcadim = 1:30
	computePCs = 50
	nneighbors = 30
	varGeneNum = 3000
}

library(reshape2)
library(Seurat)
library(dplyr)
library(stringr)

print("Start")
setwd(outpath)
#根据输入文件不同，会有不同的读数据方法

allcounts = Read10X(data.dir=countsFile,gene.column = genecolumn)

dir.create("1.Normalization")
dir.create("2.PCAAnalysis")

#将counts表中表达量为0的列删掉
allcounts = allcounts[,colSums(allcounts)>0,drop=F]

mito.genes = read.table(paste0("/MT.txt"), sep = "\t", header = FALSE)[,1]
mito.genes = mito.genes[is.finite(match(mito.genes,rownames(allcounts)))]
percent.mito <- Matrix::colSums(allcounts[mito.genes, ])/Matrix::colSums(allcounts)

#构建seuratobject
seuset <- CreateSeuratObject(counts = allcounts,project = prefix,min.cells = geneExpMinCell, min.features = 1)

print("Remove MTGene")
seuset = subset(seuset,features = setdiff(row.names(seuset),mito.genes))
seuset <- AddMetaData(object = seuset, metadata = percent.mito, col.name = "percent.mito")

#添加样本信息
sampleInfo = read.table(samplelist, sep = "\t", header = FALSE,row.names=1)
samples = sapply(strsplit(colnames(seuset), "-"), function(v) return(v[2]))
sampleInfo = sampleInfo[samples,]

names(sampleInfo) = colnames(seuset)
seuset <- AddMetaData(object = seuset, metadata = sampleInfo, col.name = "orig.ident")
write.table(data.frame(CellName=rownames(seuset@meta.data),seuset@meta.data),file="1.Normalization/CellInfosRaw.txt",sep="\t",quote=F,row.names=F)

print(paste0("nGeneFilterRegion:",cellMinGene,"-",cellMaxGene))
seuset <- subset(x = seuset, subset = nFeature_RNA > cellMinGene & nFeature_RNA < cellMaxGene)
seuset <- subset(x = seuset, subset = percent.mito < MaxMTPercent)

seuset <- FindVariableFeatures(seuset, selection.method = "vst", nfeatures = varGeneNum)

ccRDS = "./CCGenes.rds";
ccgene = readRDS(ccRDS)
seuset = CellCycleScoring(seuset, s.features =ccgene$s.genes, g2m.features = ccgene$g2m.genes,  set.ident = FALSE)
seuset$CC.Difference <- seuset$S.Score - seuset$G2M.Score

#标化为norm&scale时，scaledata
if(NormType=="scale"){
	top10 <- head(VariableFeatures(seuset), 10)
	write.table(data.frame(Gene=rownames(seuset@hvg.info),seuset@hvg.info),file="1.Normalization/hvgInfo.txt",sep="\t",row.names=F,quote=F)
	write.table(VariableFeatures(seuset),file="1.Normalization/varGene.txt",sep="\t",row.names=F,quote=F)
   	scalegene = VariableFeatures(seuset)
	print(paste0("Use Variable genes to scale data"))
	seuset <- ScaleData(object = seuset, features = scalegene, vars.to.regress = regressOut, model.use = ScaleModel)
}
#标化类型为sctransform，运行SCTransform，assay为SCT，
if(NormType=="sctransform"){
	seuset <- SCTransform(seuset, vars.to.regress = regressOut, variable.features.n = varGeneNum, return.only.var.genes=RunFast)
	top10 <- head(VariableFeatures(seuset), 10)
	write.table(VariableFeatures(seuset),file="1.Normalization/varGene.txt",sep="\t",row.names=F,quote=F)
}
write.table(data.frame(CellName=rownames(seuset@meta.data),seuset@meta.data),file="1.Normalization/CellInfosFilter.txt",sep="\t",quote=F,row.names=F)

seuset <- RunPCA(object = seuset, features = VariableFeatures(seuset), npcs = computePCs, do.print = TRUE, ndims.print = 1:computePCs, nfeatures.print = 10,seed.use=2021)

pcgene <- c()
for(i in 1:computePCs){
    topgenes = TopFeatures(object = seuset[["pca"]], dim = i,balanced=TRUE, nfeatures=60)
    topgenes = data.frame(PC=paste0("PC",i),topgenes$positive,topgenes$negative)
    pcgene = rbind(pcgene,topgenes)
}
write.table(pcgene,file="2.PCAAnalysis/pcGene.txt",sep="\t",quote=F,row.names=F,col.names=T)
write.table(data.frame(CellName=rownames(seuset@reductions$pca@cell.embeddings),seuset@reductions$pca@cell.embeddings),file="2.PCAAnalysis/pca.txt",sep="\t",quote=F,row.names=F)
print("RunUMAP")
seuset <- RunUMAP(seuset, dims = pcadim, n.neighbors = nneighbors, reduction.name = "umap3d", n.components = 3, reduction.key = "umap3d_")
seuset <- RunUMAP(seuset, dims = pcadim, n.neighbors = nneighbors)
print("RunTSNE")
seuset <- RunTSNE(object = seuset,dims = pcadim,reduction.name = "tsne3d",dim.embed=3, reduction.key = "tSNE3d_",check_duplicates = FALSE)
seuset <- RunTSNE(object = seuset,dims = pcadim,check_duplicates = FALSE)
#保存RDS	
saveRDS(seuset, file = paste0("../",prefix,".seuset.rds"))
