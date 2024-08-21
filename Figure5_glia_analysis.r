library(URD)
library(Seurat)

Idents(drg)<-'CellType'

drg <- subset(drg1,idents=c('glia cell progenitor','satellite glia cell','schwann cell'))

drg <- FindVariableFeatures(drg,nfeatures = 2000)

var_gene <- VariableFeatures(drg)

knitr::opts_chunk$set(echo = TRUE)
rgl::setupKnitr()

obj.urd <- new("URD")

obj.urd@logupx.data = drg1@assays$RNA@data
obj.urd@meta <- drg1@meta.data[,c('nCount_RNA','nFeature_RNA','percent.mt')]
obj.urd@group.ids <- drg1@meta.data[,c('GW','CellType3','SID','louvain_2.00')]

obj.urd@group.ids <- drg1@meta.data[,c('GW','CellType3','SID','louvain_2.00','CellType5')]

obj.urd@var.genes <- drg1@assays$RNA@var.features

obj.urd@tsne.y <- as.data.frame(drg1@reductions$umap@cell.embeddings)
colnames(obj.urd@tsne.y) = c("tSNE1", "tSNE2")


obj.urd@pca.load <- as.data.frame(Loadings(drg1, reduction = "scanorama"))
obj.urd@pca.scores <- as.data.frame(drg1@reductions$scanorama@cell.embeddings)
obj.urd@pca.sdev <- as.numeric(apply(drg1@reductions$scanorama@cell.embeddings, 2, stats::sd))
obj.urd@pca.sig <- c(rep(TRUE, 60))

obj.urd <- calcDM(obj.urd, knn=100,genes.use =var_gene)

plotDimArray(obj.urd,reduction.use = 'dm',dims.to.plot = 1:2,label='louvain_2.00')

root.cells <- cellsInCluster(obj.urd, "louvain_2.00", "34")
# Then we run 'flood' simulations
obj.flood <- floodPseudotime(obj.urd, root.cells = root.cells, verbose=T)


obj.urd <- floodPseudotimeProcess(obj.urd, obj.flood, floods.name="pseudotime")


library(dplyr)
library(stringr)
library(cowplot)
library(grid)

tipClusters = c('0','1')

Tips = list()
tipCells = c()
for (tc in tipClusters) {
curTipCells = cellsInCluster(obj.urd, "louvain_2.00", tc)
tipCells = c(tipCells, curTipCells)
Tips[[tc]] <- curTipCells
}

obj.urd@group.ids[tipCells, "tip.clusters"] <- obj.urd@group.ids[tipCells, "louvain_2.00"]


tipName = names(table(as.factor(obj.urd@group.ids$tip.clusters)))

obj.urd@group.ids$tip.clusters <- as.factor(as.numeric(as.factor(obj.urd@group.ids$tip.clusters)))

tipNum = names(table(obj.urd@group.ids$tip.clusters))

diffusion.logistic <- pseudotimeDetermineLogistic(obj.urd, "pseudotime", optimal.cells.forward=80,
	max.cells.back=40, pseudotime.direction="<", do.plot=T, print.values=T)


biased.tm <- pseudotimeWeightTransitionMatrix(obj.urd, pseudotime='pseudotime', 
	logistic.params=diffusion.logistic, pseudotime.direction="<")


obj.urd <- urdSubset(obj.urd, cells.keep=rownames(biased.tm))


walks = simulateRandomWalksFromTips(obj.urd, tip.group.id="tip.clusters", root.cells=root.cells, 
	transition.matrix = biased.tm, n.per.tip=25000, root.visits=1, 
	verbose=T, max.steps=5000)

obj.urd <- processRandomWalksFromTips(obj.urd, walks, verbose = T)



myTree_unprocessed <- loadTipCells(obj.urd, "tip.clusters")

divergence.methods <- c("preference")


myTree = myTree_unprocessed

myTree <- buildTree(myTree,
                        pseudotime = "pseudotime",
                        tips.use=tipNum,
                        divergence.method = 'preference',
                        visit.threshold = 0.7,
                        minimum.visits = 2,
                        cells.per.pseudotime.bin = 10,
                        bins.per.pseudotime.window = 5,
                        save.all.breakpoint.info = T,
                        save.breakpoint.plots = "breakpoint_plot",
                        min.cells.per.segment = 0.01,  # 1
                        p.thresh = 0.05,
                       weighted.fusion=T)

myTree <- nameSegments(myTree, segments=tipNum, segment.names = c('0','1'), short.names = c('0','1'))




