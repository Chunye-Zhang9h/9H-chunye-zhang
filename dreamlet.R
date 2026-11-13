library(dreamlet)
library(Seurat)
library(SingleCellExperiment)
library(ggplot2)

rdsFile = "./seuset.rds";
sampleFile = "./sample.txt";
clusterFile = "./cluster.txt";
compareFile = "./ToCompare.txt";

method = "fixed";
min_cells = 10;
min_count = 10;
min_samples = 4;
FDR = 0.05;
abs_logFC = 0.25;

# sample/cluster 文件：带表头的两列，分别为细胞 ID、样本 ID/cluster。
# compare 文件：带表头的 Case、Ctrl 两列；logFC 方向为 Case-Ctrl。
seuset = readRDS(rdsFile)
samples = read.delim(sampleFile, colClasses = "character", check.names = FALSE)
clusters = read.delim(clusterFile, colClasses = "character", check.names = FALSE)
comparisons = read.delim(compareFile, colClasses = "character", check.names = FALSE)
names(samples) = c("Cell", "Sample")
names(clusters) = c("Cell", "Cluster")

if(packageVersion("SeuratObject") >= "5.0.0"){
    counts = GetAssayData(seuset, assay = "RNA", layer = "counts")
}else{
    counts = GetAssayData(seuset, assay = "RNA", slot = "counts")
}
counts = counts[, samples$Cell]
labels = data.frame(sample_id = samples$Sample,
                    cluster_id = clusters$Cluster[match(samples$Cell, clusters$Cell)],
                    row.names = samples$Cell)
sce = SingleCellExperiment(assays = list(counts = counts), colData = S4Vectors::DataFrame(labels))
sce$sample_id = factor(sce$sample_id)
sce$cluster_id = factor(sce$cluster_id)
pb = aggregateToPseudoBulk(sce, assay = "counts", sample_id = "sample_id", cluster_id = "cluster_id")
cell_counts = cellCounts(pb)
for(i in seq_len(nrow(comparisons))){
    compare_pair = c(comparisons$Case[i], comparisons$Ctrl[i])

    # 配对模型仅保留两个 cluster 均达到细胞数阈值的样本。
    included = cell_counts[colnames(pb), compare_pair, drop = FALSE] >= min_cells
    if(method == "none"){
        selected = rowSums(included) > 0
    }else{
        selected = rowSums(included) == 2
    }
    fit = dreamletCompareClusters(pb[, selected, drop = FALSE], assays = compare_pair,
                                  method = method, collapse = TRUE, min.cells = min_cells,
                                  min.count = min_count, min.samples = min_samples,
                                  normalize.method = "TMM", BPPARAM = BiocParallel::SerialParam())

    allgenes = as.data.frame(topTable(fit, coef = "compare", number = Inf, adjust.method = "BH"))
    allgenes = cbind(Gene = rownames(allgenes), allgenes, row.names = NULL)
    allgenes$upORdown = ifelse(allgenes$logFC > 0, "up", ifelse(allgenes$logFC < 0, "down", "unchanged"))
    significant = !is.na(allgenes$adj.P.Val) & is.finite(allgenes$logFC) &
                  allgenes$adj.P.Val <= FDR & abs(allgenes$logFC) >= abs_logFC
    filteredgenes = allgenes[significant, ]


    volcano = plotVolcano(fit, coef = "compare", cutoff = FDR)
    print(volcano)
}
