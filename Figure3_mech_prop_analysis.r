library(Seurat)

library(reticulate)

scanorama <- import('scanorama')


load(file='/DRG/adata.all_cell_Scanorama.60.rds')

#subclustering of mechanoreceptor and proprioceptors

#DimPlot(drg,group.by = 'louvain_2.00',label=T)

Idents(drg)<-'CellType'

mech <- subset(drg,idents='mechanoreceptor')

prop <- subset(drg,idents='proprioceptor')


#batch correction for mechanoreceptor

mech$batch = mech$GW

table(mech$batch)

mech.list <- SplitObject(mech,split.by = 'batch')

mech.list <- lapply(X = mech.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = mech.list,nfeatures = 2000)


length(features)

seuratobjetclist <- SplitObject(mech,split = 'batch')

assaylist <- list()
genelist <- list()
for(i in 1:length(seuratobjetclist))
{
   assaylist[[i]] <- t(as.matrix(GetAssayData(seuratobjetclist[[i]], "data"))[features,])
   genelist[[i]] <- features
}


integrated.corrected.data <- scanorama$correct(assaylist, genelist, return_dimred=TRUE, return_dense=TRUE,sigma=100,knn=as.integer(50),
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

pan.seurat[['RNA']] <- CreateAssayObject(data  = as.matrix(mech@assays$RNA@data))


DefaultAssay(pan.seurat)<-'scanorama'

pan.seurat <- FindNeighbors(pan.seurat, dims = 1:8) 
FindClusters(pan.seurat)
#FindAllMarkers(pan.seurat)



pan.seurat <- RunUMAP(pan.seurat,dims= c(1:8),min.dist = 0.22,n.components = 3)

DefaultAssay(pan.seurat) <-'RNA'


#A delta: PLAGL1,PRNP,KCNAB1,PRKAR1A,MEF2C,SCN7A,FST,PDE1C，DGKZ，CPNE8,CXCL12,FLT1
#A beta filed:S100B,SFRP1,THY1，KCNIP3，APLP2，CPNE6，NEFL，ATP1B1
#A beta field: RET,POU6F2，DLG2，COL25A1,SORL1，LGI2,EPHA4,NTRK1,TAFA1,CHRM2,EPHA4,OMG
#A beta RA:SGK1,BAIAP2L1,VSNL1,UST,MEF2C,SGPP2,SGK1,CALN1,CACNG5
#A beta RA:CALN1,CACNG5,MEF2C,UST,KCNA1

DefaultAssay(pan.seurat)<-'RNA'

FeaturePlot(pan.seurat,features = 'FLT1',order=T)

DimPlot(pan.seurat,group.by = 'GW')

mech_h <- pan.seurat

mech_h <- FindNeighbors(mech_h,dims = 1:20)

DefaultAssay(mech_h)<-'scanorama'

mech_h <- FindClusters(mech_h,resolution = 0.5)



#batch correction for proprioceptor
table(prop$GW)

prop$batch = prop$GW

table(prop$batch)

prop.list <- SplitObject(prop,split.by = 'batch')

prop.list <- lapply(X = prop.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = prop.list,nfeatures = 2000)


length(features)

seuratobjetclist <- SplitObject(prop,split = 'batch')

assaylist <- list()
genelist <- list()
for(i in 1:length(seuratobjetclist))
{
   assaylist[[i]] <- t(as.matrix(GetAssayData(seuratobjetclist[[i]], "data"))[features,])
   genelist[[i]] <- features
}

integrated.corrected.data <- scanorama$correct(assaylist, genelist, return_dimred=TRUE, return_dense=TRUE,sigma=300,knn=as.integer(100),
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

pan.seurat[['RNA']] <- CreateAssayObject(data  = as.matrix(prop@assays$RNA@data))


DefaultAssay(pan.seurat)<-'scanorama'

pan.seurat <- FindNeighbors(pan.seurat, dims = 1:10) 
FindClusters(pan.seurat)
#FindAllMarkers(pan.seurat)


p1:NTRK3,SOX5
P2:CRTAC1,TCERG1L,NFIA
P3:ETV1,PVALB,TACR3,ITGA2

pan.seurat <- RunUMAP(pan.seurat,dims= 1:10,min.dist = 0.35)
prop_h = pan.seurat




#DEGs in ion channels 
ion_channel = c('ASIC1A','ASIC1B','ASIC2A','ASIC2B','ASIC3','ASIC4',
               'PIEZO1','PIEZO2','KCNK2','KCNK4','KCNK10',
               'CACNA1H','KCNQ4','SCN11A','TRPV1','TRPV2','TRPV4',
               'TRPC1','TRPC3','TRPC5','TRPC6',
               'TRPP1','TRPP2','TRPP3','TRPA1',
               'TRPM3','TRPM4','TRPM7','TRPY1',
               'SLP3','STOM','TRPN1','MEC-4','MEC-10','MEC-2',
               'SCN1A','KCNC4','KCNA2','KCND1','KCNA1','KCNB1',
                'KCNH2','KCNQ2','KCNC1','CACNA1H',
               'SCN1A','SCN8A','KCNA1','KCNA2','KCNH2','KCND1',
                'KCNC4','KCNQ2','HCN1','CACNA1A',
               'SCN1A','SCN8A','KCNA1','KCNA2','KCNH2','KCND1','KCNQ2',
                'KCNB2','KCNC4','HCN1','SCN1A','KCNA1','KCND1','KCNA2',
                'KCNH2','KCNQ2','KCNC4','KCNB1','CACNA1A','HCN1',
               'SCN1A','SCN8A','KCNC1','KCNH2','KCNA2','KCNC4','KCNA1',
                'KCNC3','KCNB2','KCND1')

#selected
ion_channel = strsplit(ion_channel,',')

ion_channel = do.call(rbind,ion_channel)

ion_channel = as.data.frame(ion_channel)

ion_channel = (unique(ion_channel$V1))

#k-channel:

library(ggplot2)
mech$set = 'mech'
prop$set = 'prop'
mech_prop = merge(mech_h,prop_h)

DefaultAssay(mech_prop)
Idents(mech_prop) = 'set
mech_prop_deg <- FindAllMarkers(mech_prop,only.pos = T)



#select ion channels in mech-prop DEGs
pos <- which(c(c('ASIC1','ASIC2','CACNG5','CACNA1A','CACHD1','CACNA1D','CACNA2D1',
                                   
                                   'CACNA1C','CACNA1A','KCNT2','KCNMB2','KCNH7','KCTD6',
                                    'KCNV1','KCNK2','KCNK4','KCNQ4','KCNN2','KCNJ3','KCNIP4','KCNU1',
                                    'KCTD3','KCNJ6','KCNQ2','KCNC1','KCNN1','KCNB1','KCNH5','KCNH8','KCTD8',
                                    'KCNB2','KCND3','KCND2','KCNH1','KCNQ5','SCN10A','SCN5A','SCN11A','SCN7A','SCN2A',
                                    'SCN8A','SCN3B','SCN1A','SCN11A','TRPC5','TRPC1','TRPC3','TRPC5','TRPM4','TRPM3','TRPV1',
                                    'PIEZO1','PIEZO2'))%in%mech_prop_deg$gene)



#mechanoreceptor URD analysis

library(URD)

DefaultAssay(mech_h)<-'scanorama'

mech_h <- FindClusters(mech_h,resolution = 0.5)

mech_h$CellType = mech_h$mech_type

DefaultAssay(mech_h)<-'scanorama'

mech_h <- FindVariableFeatures(mech_h,nfeatures = 2000)
length(VariableFeatures(mech_h))

knitr::opts_chunk$set(echo = TRUE)
rgl::setupKnitr()

drg <- subset(drg,cells = colnames(mech_h))

all(colnames(drg)==colnames(mech_h))

drg$mech_type = mech_h$mech_type
drg$scanorama_snn_res.5 = mech_h$scanorama_snn_res.5

mech_h = drg

obj.urd <- new("URD")

obj.urd@logupx.data = mech_h@assays$RNA@data
obj.urd@meta <- mech_h@meta.data[,c('nCount_RNA','nFeature_RNA','percent.mt')]
obj.urd@group.ids <- mech_h@meta.data[,c('GW','CellType3','SID','RNA_snn_res.1','mech_type')]



mech_h <- FindVariableFeatures(mech_h)

obj.urd@var.genes <- mech_prop@assays$RNA@var.features

obj.urd@var.genes <- mech_h@assays$RNA@var.features

obj.urd@tsne.y <- as.data.frame(mech_h@reductions$umap@cell.embeddings)
colnames(obj.urd@tsne.y) = c("tSNE1", "tSNE2")


Loadings(mech_h,reduction = 'scanorama')

obj.urd@pca.load <- as.data.frame(Loadings(mech_h, reduction = 'scanorama'))
obj.urd@pca.scores <- as.data.frame(mech_h@reductions$scanorama@cell.embeddings)
obj.urd@pca.sdev <- as.numeric(apply(mech_h@reductions$scanorama@cell.embeddings, 2, stats::sd))
obj.urd@pca.sig <- c(rep(TRUE, 50))


var_gene <- c(VariableFeatures(mech_h))
length(var_gene)

obj.urd <- calcDM(obj.urd, knn=100, sigma.use=30)

mech_h <- RunPCA(mech_h)

mech_h <- FindNeighbors(mech_h)

mech_h <- FindClusters(mech_h,resolution = 1)



DimPlot(mech_h,group.by = 'RNA_snn_res.1',label=T)

plotDimArray(obj.urd,reduction.use = 'dm',dims.to.plot = 1:2,label='RNA_snn_res.1')

plotDim(obj.urd, "mech_type", transitions.plot = 10000, plot.title="Clusters with transitions")

#correspond to precursor cells
root.cells <- cellsInCluster(obj.urd, "RNA_snn_res.1", "9")
# Then we run 'flood' simulations
obj.flood <- floodPseudotime(obj.urd, root.cells = root.cells, verbose=T)
# The we process the simulations into a pseudotime


obj.urd <- floodPseudotimeProcess(obj.urd, obj.flood, floods.name="pseudotime")


pseudotimePlotStabilityOverall(obj.urd)


plotDists(obj.urd, "pseudotime", "GW" ,plot.title="Pseudotime by stage")

library(dplyr)
library(stringr)
library(cowplot)
library(grid)


tipClusters = c('M_2','M_3','M_4')

Tips = list()
tipCells = c()
for (tc in tipClusters) {
curTipCells = cellsInCluster(obj.urd, "mech_type", tc)
tipCells = c(tipCells, curTipCells)
Tips[[tc]] <- curTipCells
}

obj.urd@group.ids$tip.clusters=NULL

obj.urd@group.ids[tipCells, "tip.clusters"] <- obj.urd@group.ids[tipCells, "mech_type"]


tipName = names(table(as.factor(obj.urd@group.ids$tip.clusters)))

obj.urd@group.ids$tip.clusters <- as.factor(as.numeric(as.factor(obj.urd@group.ids$tip.clusters)))

tipNum = names(table(obj.urd@group.ids$tip.clusters))

print(tipNum)
print(tipName)

diffusion.logistic <- pseudotimeDetermineLogistic(obj.urd, "pseudotime", optimal.cells.forward=20,
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



cells.per.pseudotime.bins <- c(80)  # 80
bins.per.pseudotime.windows <- c(5)  # 5
visit.thresholds <- c(0.7)  # 0.7
minimum.visitss <- c(10)  # 10


myTree = myTree_unprocessed

myTree <- buildTree(myTree,
                        pseudotime = "pseudotime",
                        tips.use=tipNum,
                        divergence.method = 'preference',
                        visit.threshold = 0.05,
                        minimum.visits = 10,
                        cells.per.pseudotime.bin = 65,
                        bins.per.pseudotime.window = 5,
                        save.all.breakpoint.info = T,
                        save.breakpoint.plots = "breakpoint_plot",
                        min.cells.per.segment = 0.01,  # 1
                        p.thresh = 0.05,
                       weighted.fusion=T)

myTree <- nameSegments(myTree, segments=tipNum, segment.names = tipNum, short.names = tipNum)
                 

#proprioceptor URD analysis

library(URD)


knitr::opts_chunk$set(echo = TRUE)
rgl::setupKnitr()

drg <- subset(drg,cells = colnames(prop_h)))

all(colnames(drg)==colnames(prop_h))

colnames(prop_h@meta.data)

DimPlot(prop_h,group.by = 'scanorama_snn_res.1',label=T)

drg$prop_type = prop_h$prop_type
drg$scanorama_snn_res.1 = prop_h$scanorama_snn_res.1

prop_h = drg1

obj.urd <- new("URD")

obj.urd@logupx.data = prop_h@assays$RNA@data
obj.urd@meta <- prop_h@meta.data[,c('nCount_RNA','nFeature_RNA','percent.mt')]
obj.urd@group.ids <- prop_h@meta.data[,c('GW','CellType3','SID','scanorama_snn_res.1','prop_type')]



prop_h <- FindVariableFeatures(prop_h,nfeatures=1500)

obj.urd@var.genes <- prop_h@assays$RNA@var.features

obj.urd@tsne.y <- as.data.frame(prop_h@reductions$umap@cell.embeddings)
colnames(obj.urd@tsne.y) = c("tSNE1", "tSNE2")


obj.urd@pca.load <- as.data.frame(Loadings(prop_h, reduction = 'scanorama'))
obj.urd@pca.scores <- as.data.frame(prop_h@reductions$scanorama@cell.embeddings)
obj.urd@pca.sdev <- as.numeric(apply(prop_h@reductions$scanorama@cell.embeddings, 2, stats::sd))
obj.urd@pca.sig <- c(rep(TRUE, 50))

var_gene <- c(VariableFeatures(prop_h))
length(var_gene)

obj.urd <- calcDM(obj.urd, knn=100, sigma.use=30)

plotDimArray(obj.urd,reduction.use = 'dm',dims.to.plot = 1:2,label='scanorama_snn_res.1')

plotDim(obj.urd, "prop_type", transitions.plot = 10000, plot.title="Clusters with transitions")

#corresponding to precursor cells
root.cells <- cellsInCluster(obj.urd, "scanorama_snn_res.1", "12")
# Then we run 'flood' simulations
obj.flood <- floodPseudotime(obj.urd, root.cells = root.cells, verbose=T)
# The we process the simulations into a pseudotime


obj.urd <- floodPseudotimeProcess(obj.urd, obj.flood, floods.name="pseudotime")


pseudotimePlotStabilityOverall(obj.urd)

plotDists(obj.urd, "pseudotime", "GW" ,plot.title="Pseudotime by stage")

library(dplyr)
library(stringr)
library(cowplot)
library(grid)

tipClusters = c('P_2','P_3')

Tips = list()
tipCells = c()
for (tc in tipClusters) {
curTipCells = cellsInCluster(obj.urd, "prop_type", tc)
tipCells = c(tipCells, curTipCells)
Tips[[tc]] <- curTipCells
}

obj.urd@group.ids$tip.clusters=NULL

obj.urd@group.ids[tipCells, "tip.clusters"] <- obj.urd@group.ids[tipCells, "prop_type"]


tipName = names(table(as.factor(obj.urd@group.ids$tip.clusters)))

obj.urd@group.ids$tip.clusters <- as.factor(as.numeric(as.factor(obj.urd@group.ids$tip.clusters)))

tipNum = names(table(obj.urd@group.ids$tip.clusters))

print(tipNum)
print(tipName)

diffusion.logistic <- pseudotimeDetermineLogistic(obj.urd, "pseudotime", optimal.cells.forward=20,
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

tipNum

cells.per.pseudotime.bins <- c(80)  # 80
bins.per.pseudotime.windows <- c(5)  # 5
visit.thresholds <- c(0.7)  # 0.7
minimum.visitss <- c(10)  # 10


myTree = myTree_unprocessed

myTree <- buildTree(myTree,
                        pseudotime = "pseudotime",
                        tips.use=tipNum,
                        divergence.method = 'preference',
                        visit.threshold = 0.05,
                        minimum.visits = 10,
                        cells.per.pseudotime.bin = 65,
                        bins.per.pseudotime.window = 5,
                        save.all.breakpoint.info = T,
                        save.breakpoint.plots = "breakpoint_plot",
                        min.cells.per.segment = 0.01,  # 1
                        p.thresh = 0.05,
                       weighted.fusion=T)



