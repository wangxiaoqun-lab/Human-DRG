library(Seurat)

library(reticulate)

#subset nociceptors

#subclustering of nociceptors
noci = subset(drg,idents='nociceptor')
noci$batch = noci$GW_Rep
seuratobjetclist = unique(noci$batch)

seuratobjetclist <- SplitObject(noci.combined,split = 'batch')
noci = FindVariableFeatures(noci,nfeatures=2000)
var_gene = VariableFeatures(noci)
assaylist <- list()
genelist <- list()
for(i in 1:length(seuratobjetclist))
{
   assaylist[[i]] <- t(as.matrix(GetAssayData(seuratobjetclist[[i]], "data"))[var_gene,])
   genelist[[i]] <- var_gene
}

#integrated.data <- scanorama$integrate(assaylist, genelist,return_dimred=TRUE,return_dense=TRUE)
#corrected.data <- scanorama$correct(assaylist, genelist, return_dense=TRUE)
integrated.corrected.data <- scanorama$correct(assaylist, genelist, return_dimred=TRUE, return_dense=TRUE,sigma=100,knn=as.integer(200),
                                              dimred=as.integer(50))

intdata <- lapply(integrated.corrected.data[[2]], t)
panorama <- do.call(cbind, intdata)


rownames(panorama) <- as.character(integrated.corrected.data[[3]])
colnames(panorama) <- unlist(sapply(assaylist, rownames))

intdimred <- do.call(rbind, integrated.corrected.data[[1]])
colnames(intdimred) <- paste0("PC_", 1:50)

stdevs <- apply(intdimred, MARGIN = 2, FUN = sd)

pan.seurat <- CreateSeuratObject(counts = panorama, assay = "scanorama",  project = "human_noci")

pan.seurat@meta.data <- do.call(rbind, lapply(seuratobjetclist, function(x) x@meta.data))

rownames(pan.seurat@meta.data) <- colnames(pan.seurat)
rownames(intdimred) <- colnames(pan.seurat)


pan.seurat[["pca"]] <- CreateDimReducObject(embeddings = intdimred, stdev = stdevs, key = "PC_", assay = "scanorama")

pan.seurat[['RNA']] <- CreateAssayObject(data  = as.matrix(noci.combined@assays$RNA@data))


DefaultAssay(pan.seurat)<-'RNA'

table(pan.seurat$pp_np_type4)

pan.seurat <- FindNeighbors(pan.seurat, dims = 1:8) 
FindClusters(pan.seurat)
#FindAllMarkers(pan.seurat)

pan.seurat <- RunUMAP(pan.seurat,dims= 1:8,min.dist = 0.3)

DefaultAssay(pan.seurat)<-'scanorama'

pan.seurat <- FindClusters(pan.seurat,resolution = 2)






#construction of nociceptor URD tree
library(URD)

library(harmony)

noci <- FindVariableFeatures(noci,nfeatures = 2500)

knitr::opts_chunk$set(echo = TRUE)
rgl::setupKnitr()

obj.urd <- new("URD")

obj.urd@logupx.data = noci@assays$RNA@data
obj.urd@meta <- noci@meta.data[,c('nCount_RNA','nFeature_RNA','percent.mt')]
obj.urd@group.ids <- noci@meta.data[,c('GW','CellType3','SID','louvain_2.00','subtype')]

obj.urd@var.genes <- noci@assays$RNA@var.features

obj.urd@tsne.y <- as.data.frame(noci@reductions$umap@cell.embeddings)
colnames(obj.urd@tsne.y) = c("tSNE1", "tSNE2")


Loadings(noci,reduction = 'scanorama')

obj.urd@pca.load <- as.data.frame(Loadings(noci, reduction = "scanorama"))
obj.urd@pca.scores <- as.data.frame(noci@reductions$scanorama@cell.embeddings)
obj.urd@pca.sdev <- as.numeric(apply(noci@reductions$scanorama@cell.embeddings, 2, stats::sd))
obj.urd@pca.sig <- c(rep(TRUE, 60))


obj.urd <- calcKNN(obj.urd, nn=100)


var_gene <- VariableFeatures(noci)

obj.urd <- calcDM(obj.urd, knn=100, sigma.use=30,genes.use = var_gene)

plotDimArray(obj.urd,reduction.use = 'dm',dims.to.plot = c(2,3),label='subtype')

plotDim(obj.urd, "louvain_2.00", transitions.plot = 10000, plot.title="Clusters with transitions")

root.cells <- cellsInCluster(obj.urd, "louvain_2.00", "2")
# Then we run 'flood' simulations
obj.flood <- floodPseudotime(obj.urd, root.cells = root.cells, verbose=T)
# The we process the simulations into a pseudotime


obj.urd <- floodPseudotimeProcess(obj.urd, obj.flood, floods.name="pseudotime")

pseudotimePlotStabilityOverall(obj.urd)

library(dplyr)
library(stringr)
library(cowplot)
library(grid)

tipClusters = c('C-LTMR','NP2(ERBB4_OPRM1_CALCA)','NP1(GFRA2_MRGPRD)','PEP4(TRPM8)',
               'PEP5(TAFA1_NEFH_KIT)','PEP1(TAC1)','PEP3(TAC1_NTRK2)','PEP2(TAFA1_NEFH)',
               'NP3(NPPB_SST)', 'NP4(GFRA1_LPAR1)','ERBB4_CHRM2_IL33_NOS1',
               'TAC1_NOS1','GFRA2_2')

Tips = list()
tipCells = c()
for (tc in tipClusters) {
curTipCells = cellsInCluster(obj.urd, "subtype", tc)
tipCells = c(tipCells, curTipCells)
Tips[[tc]] <- curTipCells
}

obj.urd@group.ids$tip.clusters=NULL

obj.urd@group.ids[tipCells, "tip.clusters"] <- obj.urd@group.ids[tipCells, "subtype"]


tipName = names(table(as.factor(obj.urd@group.ids$tip.clusters)))

obj.urd@group.ids$tip.clusters <- as.factor(as.numeric(as.factor(obj.urd@group.ids$tip.clusters)))

tipNum = names(table(obj.urd@group.ids$tip.clusters))

diffusion.logistic <- pseudotimeDetermineLogistic(obj.urd, "pseudotime", optimal.cells.forward=40,
	max.cells.back=80, pseudotime.direction="<", do.plot=T, print.values=T)


biased.tm <- pseudotimeWeightTransitionMatrix(obj.urd, pseudotime='pseudotime', 
	logistic.params=diffusion.logistic, pseudotime.direction="<")


obj.urd <- urdSubset(obj.urd, cells.keep=rownames(biased.tm))

walks = simulateRandomWalksFromTips(obj.urd, tip.group.id="tip.clusters", root.cells=root.cells, 
	transition.matrix = biased.tm, n.per.tip=25000, root.visits=1, 
	verbose=T, max.steps=5000)


obj.urd <- processRandomWalksFromTips(obj.urd, walks, verbose = T)

myTree_unprocessed <- loadTipCells(obj.urd, "tip.clusters")


myTree <- buildTree(myTree,
                    pseudotime = "pseudotime",
                    tips.use=tipNum,
                    divergence.method = 'preference',
                    visit.threshold = 0.05,
                    minimum.visits = 5,
                    cells.per.pseudotime.bin =30,
                    bins.per.pseudotime.window = 1,
                    save.all.breakpoint.info = T,
                    save.breakpoint.plots = "breakpoint_plot",
                    min.cells.per.segment =50,  # 1
                    p.thresh = 0.05,
                    weighted.fusion=T,
                    min.pseudotime.per.segment = .05,
                    dendro.cell.dist.to.tree=0.1)

myTree <- nameSegments(myTree, segments=tipNum, segment.names = tipNum, short.names = tipNum)
                
plotTree(myTree,'GW')+scale_color_manual(values=c('GW7'='#39294f','GW8'='#4d588e',
                                                         'GW9'='#5f6f92','GW10'='#739b96',
                                                         'GW12'='#abc6b4','GW14'='#9dbf65',
                                                         'GW15'='#c2d04e','GW17'='#f3e767',
                                                         'GW21'='yellow'))

