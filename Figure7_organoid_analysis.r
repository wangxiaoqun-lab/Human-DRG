library(Seurat)

#import d30-1 sample
d30_1 = Read10X_h5(file='/newdata/mengdi_mount/DRG/organoid_seq/filtered_gene_bc_matrices_h5.h5')

d30_1 = CreateSeuratObject(d30_1,min.cells = 3,min.features = 0,project = 'd30_1')

d30_1 = NormalizeData(d30_1)

d30_1 = FindVariableFeatures(d30_1)

d30_1 = ScaleData(d30_1)

d30_1 = RunPCA(d30_1)

d30_1 = RunUMAP(d30_1,dims=1:15)


#import d30-2 sample

d30_2 = Read10X_h5(file='/newdata/mengdi_mount/DRG/organoid_seq/d30_2_filtered_gene_bc_matrices_h5.h5')

d30_2 = CreateSeuratObject(d30_2,min.cells = 3,min.features = 0,project = 'd30_2')

d30_2 = NormalizeData(d30_2)

d30_2 = FindVariableFeatures(d30_2)

d30_2 = ScaleData(d30_2)

d30_2 = RunPCA(d30_2)

d30_2 = RunUMAP(d30_2,dims = 1:30)

#import d60-1 sample


d60_1 = Read10X_h5(file='/newdata/mengdi_mount/DRG/organoid_seq/d60_1filtered_gene_bc_matrices_h5.h5')

d60_1 = CreateSeuratObject(d60_1,min.cells = 3,min.features = 0,project = 'd60_1')

d60_1 = NormalizeData(d60_1)
d60_1 = FindVariableFeatures(d60_1)

d60_1 = ScaleData(d60_1)

d60_1 = RunPCA(d60_1)

d60_1 = RunUMAP(d60_1,dims = 1:30)

#import d60-2 sample

d60_2 = Read10X_h5(file='/newdata/mengdi_mount/DRG/organoid_seq/d60_2filtered_gene_bc_matrices_h5.h5')

d60_2 = CreateSeuratObject(d60_2,min.cells = 3,min.features = 0,project = 'd60_2')

d60_2 = NormalizeData(d60_2)
d60_2 = FindVariableFeatures(d60_2)

d60_2 = ScaleData(d60_2)
d60_2 = RunPCA(d60_2)

d60_2 = RunUMAP(d60_2,dims = 1:30)

#merge all 4 samples and performed QC 
org = merge(d30,d60_1)

org = merge(org,d60_2)

org[["percent.mt"]] <- PercentageFeatureSet(org, pattern = "^MT-")

org$set = 'organoid'

VlnPlot(org, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,group.by = 'set')

org <- subset(org, subset = nFeature_RNA > 300 & nFeature_RNA < 5000&nCount_RNA<20000 & percent.mt < 5)

org = JoinLayers(org)

org <- NormalizeData(org, normalization.method = "LogNormalize", scale.factor = 10000)

org = FindVariableFeatures(org)

org = ScaleData(org)

org = RunPCA(org)

org = RunUMAP(org,dims=1:25)

org = FindNeighbors(org,dims = 1:25)

org = FindClusters(org,resolution = 1)

#batch correction 
org = JoinLayers(org)

org[["RNA"]] <- split(org[["RNA"]], f = org$orig.ident)


org <- NormalizeData(org)
org <- FindVariableFeatures(org)
org <- ScaleData(org)
org <- RunPCA(org)

org <- FindNeighbors(org, dims = 1:30, reduction = "pca")
org <- FindClusters(org, resolution = 2, cluster.name = "unintegrated_clusters")

org <- IntegrateLayers(
  object = org, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = FALSE,k.weight=150
)

org <- FindNeighbors(org, reduction = "integrated.cca", dims = 1:30)
org <- FindClusters(org, resolution = 2, cluster.name = "cca_clusters")

library(Signac)

org <- RunUMAP(org, reduction = "integrated.cca", dims = c(1:25), reduction.name = "umap.cca1",n.components = 2,min.dist = 0.5) 

#check sequencing depth correlation
DepthCaor(org,reduction = 'integrated.cca',n=30)

#dimension reduction analysis
org <- RunUMAP(org, reduction = "integrated.cca", dims = c(2,3,5:25), reduction.name = "umap.cca1",n.components = 2,min.dist = 0.5) 

org <- FindNeighbors(org, reduction = "integrated.cca", dims = c(2,3,5:25))
org <- FindClusters(org, resolution = 2, cluster.name = "cca_clusters")



#subclustering of sensory neurons in hDRGOs

org_f_N = subset(org_f,idents='sensory neuron')

org_f_N = RunUMAP(org_f_N,dims =1:15,reduction = 'integrated.cca')




#integration of human DRG with hDRGOs

Idents(drg)='CellType'

drg = subset(drg,idents = c('muscle cell','vascular smooth muscle cell','red blood cells','ribo'),invert=T)

drg$set = 'human'

org$set = 'org'

merge = merge(drg,org)

pos = grep('human',merge$set)
merge$orig.ident[pos] = 'human'

pos = grep('human',merge$set)
merge$orig.ident[pos] = merge$GW[pos]

merge$batch = merge$orig.ident

pos = grep('GW7|GW8|GW9',merge$batch)
merge$batch[pos] = 'batch1'
pos = grep('GW10|GW12',merge$batch)
merge$batch[pos] = 'batch2'
pos = grep('GW14|GW15',merge$batch)
merge$batch[pos] = 'batch3'
pos = grep('GW17|GW21',merge$batch)
merge$batch[pos] = 'batch4'


merge = JoinLayers(merge)



merge[["RNA"]] <- split(merge[["RNA"]], f = merge$batch)

merge <- NormalizeData(merge)
merge <- FindVariableFeatures(merge)
merge <- ScaleData(merge)
merge <- RunPCA(merge)

merge <- FindNeighbors(merge, dims = 1:30, reduction = "pca")
merge <- FindClusters(merge, resolution = 2, cluster.name = "unintegrated_clusters")

merge <- RunUMAP(merge, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
# visualize by batch and cell type annotation
# cell type annotations were previously added by Azimuth


options(future.globals.maxSize = 3e+09)

merge <- IntegrateLayers(
  object = merge, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = FALSE)

merge <- FindNeighbors(merge, reduction = "integrated.cca", dims = 1:25)
merge <- FindClusters(merge, resolution = 2, cluster.name = "cca_clusters")

merge <- RunUMAP(merge, reduction = "integrated.cca", dims = 1:25, reduction.name = "umap.cca")









