library(methods)
library(velocyto.R)
library(pagoda2)

outpath = "./result";
loom = "./merge.loom";
clusterFile = "./cellSummary.txt";
coordsFile = "./coords.2d.txt"
quantile = 0.02

setwd(outpath)
ldat <- read.loom.matrices(loom)

cluster = read.table(clusterFile,sep="\t",header=T,row.names=1)

coords = read.table(coordsFile,sep="\t",header=T,row.names=1)[,1:2]

emat <- ldat$spliced
nmat <- ldat$unspliced
select = is.finite(match(colnames(emat),rownames(cluster)))
emat <- emat[,select]
nmat <- nmat[,select]

r <- Pagoda2$new(emat,modelType='plain',trim=10,log.scale=T)
r$adjustVariance(plot=F,do.par=T,gam.k=10)
r$calculatePcaReduction(nPcs=100,n.odgenes=3e3,maxit=300)
pca = r$reductions$PCA

cluster.label = cluster[,1]
names(cluster.label) = rownames(cluster)
cell.colors <- pagoda2:::fac2col(cluster.label)

emb <- as.matrix(coords)
cell.dist <- as.dist(1-armaCor(t(pca)))

emat <- filter.genes.by.cluster.expression(emat,cluster.label,min.max.cluster.average = 0.5)
nmat <- filter.genes.by.cluster.expression(nmat,cluster.label,min.max.cluster.average = 0.05)
length(intersect(rownames(emat),rownames(emat)))

rvel.cd <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells=20,cell.dist=cell.dist,fit.quantile=quantile)

pdf("Velocity.pdf",width=10,height=10)
info = show.velocity.on.embedding.cor(emb,rvel.cd,n=300,scale='sqrt',cell.colors=ac(cell.colors,alpha=0.5),cex=0.8,arrow.scale=5,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=1,do.par=F,cell.border.alpha = 0.2,return.details = TRUE)
dev.off()

