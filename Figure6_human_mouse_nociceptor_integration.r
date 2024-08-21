library(Seurat)

load(file='DRG/nature_ginty/ginty_all_neuron_clustering.Robj')
ortho_gene = read.csv(file='DRG/biomart_human_mouse_one2one_ortho.csv)

DefaultAssay(ginty)<-'RNA'

table(ginty$type)

pos <- grep('CGRP|nociceptors|SST|TRPM8',ginty$type)
table(ginty$type[pos])

mouse_noci = subset(ginty,cells = colnames(ginty)[pos])

table(mouse_noci$type)


count = as.matrix(mouse_noci@assays$RNA@counts)
rownames(count) = toupper(rownames(count))
mouse_noci1 = CreateSeuratObject(counts = count,min.cells = 0,min.features = 0,project = 'mouse_noci')

#dimension reduction
all(colnames(mouse_noci)==colnames(mouse_noci1))

mouse_noci1$time = mouse_noci$time

mouse_noci1$type = mouse_noci$type

mouse_noci1 = NormalizeData(mouse_noci1)

mouse_noci1 = FindVariableFeatures(mouse_noci1)

mouse_noci1 = ScaleData(mouse_noci1)

mouse_noci1 = RunPCA(mouse_noci1)

mouse_noci1 = FindNeighbors(mouse_noci1,dims = 1:30)

mouse_noci1 = FindClusters(mouse_noci1)

mouse_noci1 = RunUMAP(mouse_noci1,dims = 1:10)

DimPlot(mouse_noci1,group.by = 'type')

length(VariableFeatures(mouse_noci1))

var_gene <- VariableFeatures(mouse_noci1)
length(var_gene)

seuratobjetclist <- SplitObject(mouse_noci1,split = 'time')


assaylist <- list()
genelist <- list()
for(i in 1:length(seuratobjetclist))
{
   assaylist[[i]] <- t(as.matrix(GetAssayData(seuratobjetclist[[i]], "data"))[var_gene,])
   genelist[[i]] <- var_gene
}

#batch correction for Ginty's data
library(reticulate)

scanorama <- import('scanorama')

integrated.corrected.data <- scanorama$correct(assaylist, genelist, return_dimred=TRUE, return_dense=TRUE,sigma=150,knn=as.integer(100),
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

pan.seurat[['RNA']] <- CreateAssayObject(data  = (ginty@assays$RNA@data[,colnames(pan.seurat)]))


all(colnames(pan.seurat)==colnames(ginty))

DefaultAssay(pan.seurat)<-'scanorama'

pan.seurat <- FindNeighbors(pan.seurat, dims = 1:20) 
pan.seurat <- FindClusters(pan.seurat,resolution = 5)

pan.seurat_sub<- RunUMAP(pan.seurat,dims= c(1:20),min.dist = 0.3,n.components = 2)

DimPlot(pan.seurat_sub,group.by='type',label=T,dims = c(1,2),order=c('CGRP_Theta','CGRP_Gamma',
                                                                    'CGRP_Alpha','nonpeptidergic_nociceptors',
                                                                    'C-LTMR','SST','CGRP_Zeta','CGRP_Epsilion',
                                                                    'proprioceptor','Adelta_LTMR','TRPM8',
                                                                    'Abeta_RA_LTMR',
                                                                    'CGRP_Eta','Abeta_Filed'))+NoLegend()

saveRDS(pan.seurat_sub,file='/DRG/nature_ginty/ginty_all_nociceptor_scanorama_clustering.rds')


#subset data 
Idents(mouse_noci1)='time'

table(Idents(mouse_noci1))

mouse_noci2 = subset(mouse_noci1,idents=c('E15.5','P0','P5'))

Idents(mouse_noci2) = 'time'

DefaultAssay(mouse_noci2)

mouse_noci2 = FindVariableFeatures(mouse_noci2,nfeatures = 2000)

length(VariableFeatures(mouse_noci2))

var_gene <- VariableFeatures(mouse_noci2)
length(var_gene)

seuratobjetclist <- SplitObject(mouse_noci2,split = 'time')



assaylist <- list()
genelist <- list()
for(i in 1:length(seuratobjetclist))
{
   assaylist[[i]] <- t(as.matrix(GetAssayData(seuratobjetclist[[i]], "data"))[var_gene,])
   genelist[[i]] <- var_gene
}



library(reticulate)

scanorama <- import('scanorama')


integrated.corrected.data <- scanorama$correct(assaylist, genelist, return_dimred=TRUE, return_dense=TRUE,sigma=150,knn=as.integer(100),
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


count_raw = ginty@assays$RNA@data

dim(count_raw)

rownames(count_raw) = toupper(rownames(count_raw))

dim(count_raw)

pan.seurat

all(colnames(pan.seurat)%in%colnames(count_raw))

all(rownames(count_raw))

pan.seurat[["pca"]] <- CreateDimReducObject(embeddings = intdimred, stdev = stdevs, key = "PC_", assay = "scanorama")

pan.seurat[['RNA']] <- CreateAssayObject(data  = count_raw[,colnames(pan.seurat)])


all(colnames(pan.seurat)==colnames(mouse_noci2))

DefaultAssay(pan.seurat)<-'scanorama'

pan.seurat <- FindNeighbors(pan.seurat, dims = 1:20) 
pan.seurat <- FindClusters(pan.seurat,resolution = 5)

pan.seurat <- RunUMAP(pan.seurat,dims= c(1:15),min.dist = 0.3,n.components = 2)

DimPlot(pan.seurat,group.by='type',label=T)

DefaultAssay(pan.seurat) = 'RNA'

DimPlot(pan.seurat,label=T,group.by = 'type')

FeaturePlot(pan.seurat,features = c('NTRK2'),order=F)

pan.seurat
saveRDS(pan.seurat,file='/gpfs2/wulab10/DRG/nature_ginty/ginty_all_nociceptors_clustering_scanorama.rds')

mouse_ginty = pan.seurat
#subset mouse nociceptor from E12.5 to P5

mouse_noci = subset(mouse_ginty,ident=c('E12.5,'E15.5','P0','P5'))

human_noci = readRDS(file='/DRG/human_noci_subclustering.rds')

DefaultAssay(human_noci)<-'RNA'
Idents(human_noci)<-'noci_type'

#excluded C-LTMR

human_noci <- subset(human_noci,idents = 'C-LTMR',invert=T)
mouse_noci$set = 'mouse'
human_noci$set= 'human'
mouse_noci$batch = mouse_noci$time
human_noci$batch = human_noci$GW


HM_merge <- merge(mouse_noci,human_noci)
HM_merge.list <- SplitObject(HM_merge,split.by = 'batch')


#perform integration analysis
#only use one2one orthologs

features <- SelectIntegrationFeatures(HM_merge.list,nfeatures = 3000)
features = intersect(features,ortho_gene$gene)
HM_merge.anchors <- FindIntegrationAnchors(object.list = HM_merge.list,k.anchor = 5,anchor.features=features)

HM_merge.combined <- IntegrateData(anchorset = HM_merge.anchors)

DefaultAssay(HM_merge.combined)<-'integrated'

HM_merge.combined <- ScaleData(HM_merge.combined, verbose = FALSE)


DefaultAssay(HM_merge.combined)<-'integrated'

HM_merge.combined <- RunPCA(HM_merge.combined,ndims.print = 1:30)


DefaultAssay(HM_merge.combined)<-'integrated'

HM_merge.combined <- RunUMAP(HM_merge.combined,min.dist = 0.3,dims = c(1:15))

HM_merge.combined <- FindNeighbors(HM_merge.combined,dims=1:15)

HM_merge.combined <- FindClusters(HM_merge.combined,resolution = 1)

DimPlot(HM_merge.combined,label=T,group.by = 'integrated_snn_res.1',)+NoLegend()

DimPlot(HM_merge.combined,label=T,group.by = 'type')+NoLegend()

DimPlot(HM_merge.combined,group.by = 'pp_np_type6',label=T)+NoLegend()


eta:TAFA1,
Gamma:TAC1_NOS1
Alpha:TAC1
Trpm8:TRPM8
Gamma:NTRK2
Epsillion:TAFA1
zeta:DCC,ERBB4
SST:SST
non_PEP:GFRA1,GFRA2
theta:GFRA2

#quantification of cell type correspondence
Idents(HM_merge.combined)<-'set'

HM_merge.combined_m <- subset(HM_merge.combined,idents='mouse')

res_list = as.data.frame(table(HM_merge.combined$res_type1))

colnames(res_list) = c('group','number')

tmp = as.data.frame(table(HM_merge.combined$set[HM_merge.combined$res_type1=='c1'])/table(HM_merge.combined$set))

tmp = log(tmp[1,2]/tmp[2,2])

tmp = data.frame(group = 'c1',value = tmp)

tmp

for (i in (2:nrow(res_list))){
   tmp1 = as.data.frame(table(HM_merge.combined$set[HM_merge.combined$res_type1==res_list$group[i]])/table(HM_merge.combined$set)) 
   tmp1 = log(tmp1[1,2]/tmp1[2,2]) 
   tmp1 = data.frame(group = res_list$group[i],value = tmp1)
   tmp = rbind(tmp,tmp1) 
}




tmp$group = factor(tmp$group,levels = c('c1','c2','c3','c4','c5','c6','c7','c8','c9','c10','c11','c12'))

ggplot(tmp,aes(x=group,y=value))+geom_bar(stat='identity',color='black',fill='steelblue4')+geom_hline(yintercept = 2,lty='dashed',size=1,color='brown')+geom_hline(yintercept = -2,lty='dashed',size=1,color='brown')+theme_bw()+theme(panel.grid = element_blank())


pos_h = which(HM_merge.combined$set=='human')

pos_m = which(HM_merge.combined$set=='mouse')

cluster1 <- as.factor(HM_merge.combined$Subtype1[pos_h])
cluster2 <- as.factor(HM_merge.combined$Subtype1[pos_m])


names(cluster1)<-colnames(HM_merge.combined)[pos_h]

names(cluster2)<-colnames(HM_merge.combined)[pos_m]

cluster1 <- droplevels(cluster1)
  cluster2 <- droplevels(cluster2)

#HM_merge.combined$res_type = as.character(HM_merge.combined$integrated_snn_res.1)

levels(cluster_consensus)

table(Idents(HM_merge.combined))

Idents(HM_merge.combined)<-'res_type1'

levels(HM_merge.combined)<- as.factor(c('c1','c2','c3','c4','c5','c6','c7','c8','c9','c10','c11','c12'))

HM_merge.combined$res_type1 = factor(HM_merge.combined$res_type1,levels = (c('c1','c2','c3','c4','c5','c6','c7','c8','c9','c10','c11','c12')))



cluster_consensus <- as.factor(HM_merge.combined@meta.data$res_type1)


names(cluster_consensus)<- colnames(HM_merge.combined)

head(cluster_consensus)


cluster_consensus <- droplevels(cluster_consensus)

unique(cluster_consensus)

table(HM_merge.combined$res_type1)

min.frac=0.0355
min.cells=5

unique(nodes_middle)

cluster1 <- cluster1[!is.na(cluster1)]
  cluster2 <- cluster2[!is.na(cluster2)]
  nodes1 <- levels(cluster1)[table(cluster1) > 0]
  nodes2 <- levels(cluster2)[table(cluster2) > 0]
  nodes_middle <- levels(cluster_consensus)[table(cluster_consensus) > 0]
  node_Xs <- c(
    rep(1, length(nodes1)), rep(2, length(nodes_middle)),
    rep(3, length(nodes2)))

 edge_list <- list()
  for (i in 1:length(nodes1)) {
    temp <- list()
    i_cells <- names(cluster1)[cluster1 == nodes1[i]]
    for (j in 1:length(nodes_middle)) {
      if (length(which(cluster_consensus[i_cells] == nodes_middle[j])) / length(i_cells) > min.frac &
          length(which(cluster_consensus[i_cells] == nodes_middle[j])) > min.cells) {
        temp[[nodes_middle[j]]] <- sum(cluster_consensus[i_cells] ==
                                         nodes_middle[j]) / length(cluster1)
      }
    }
    edge_list[[nodes1[i]]] <- temp
  }

edge_list <- list()
  for (i in 1:length(nodes1)) {
    temp <- list()
    i_cells <- names(cluster1)[cluster1 == nodes1[i]]
    for (j in 1:length(nodes_middle)) {
      if (length(which(cluster_consensus[i_cells] == nodes_middle[j])) / length(i_cells) > min.frac &
          length(which(cluster_consensus[i_cells] == nodes_middle[j])) > min.cells) {
        temp[[nodes_middle[j]]] <- sum(cluster_consensus[i_cells] ==
                                         nodes_middle[j]) / length(cluster1)
      }
    }
    edge_list[[nodes1[i]]] <- temp
  }

cluster3 <- cluster_consensus[names(cluster2)]
  for (i in 1:length(nodes_middle)) {
    temp <- list()
    i_cells <- names(cluster3)[cluster3 == nodes_middle[i]]
    for (j in 1:length(nodes2)) {
      j_cells <- names(cluster2)[cluster2 == nodes2[j]]
      if (length(which(cluster_consensus[j_cells] == nodes_middle[i])) / length(j_cells) > min.frac &
          length(which(cluster_consensus[j_cells] == nodes_middle[i])) > min.cells) {
        if (!is.na(sum(cluster2[i_cells] == nodes2[j]))) {
          temp[[nodes2[j]]] <- sum(cluster2[i_cells] ==
                                     nodes2[j]) / length(cluster2)
        }
      }
    }
    edge_list[[nodes_middle[i]]] <- temp
  }

label.cex=0.6
label.col='black'
lab.srt=0

library(riverplot)

nodes1
nodes_middle
nodes2

'c14'='pink2',
                                                                            'c15'='#32ED3B'



node_cols <- list()
  ggplotColors <- function(g) {
    d <- 360 / g
    h <- cumsum(c(15, rep(d, g - 1)))
    grDevices::hcl(h = h, c = 100, l = 65)
  }
  pal <- c('#7FE3EA','#6856C1','#6F2E7D','#5191C0','#6A6E8D','#A79BC8',
          '#D17066','#EFE752','#A61D31','#D3489B','#ED6647','#ED356E',
           '#912D3B','#4C66AB','#7FDFE2','#CCB387','goldenrod2','#6926C2','pink2','#32ED3B',
           '#E57259','#A6940E','#833E05','#D48019','#CE6AEF','#35A489','#9B0868',
          '#B54075','royalblue','red','#D6685C','#535E7E','green','#F7C82F',
          'purple','#E5EF18','forestgreen','pink2','#2FAFF7','#7F521B')
  for (i in 1:length(nodes1)) {
    node_cols[[nodes1[i]]] <- list(col = pal[i], textcex = label.cex,
                                   textcol = label.col, srt = lab.srt)
  }
 pal <- c('#ce8600','#3c65eb','#e057ac','#e81988','tan3','yellow','lightblue1','#59be5a','#34b14e','#daf399',
          '#1f8342','#874300','#fdac34','#a79d73','red3','#d08ece','#855a83','#d08475','#a4e083','#6692b3','#9e258e','#6e7d90')
#pal <- as.character(common_pal1$pal)
  for (i in 1:length(nodes_middle)) {
    node_cols[[nodes_middle[i]]] <- list(col = pal[i], textcex = label.cex,
                                         textcol = label.col, srt = lab.srt)
  }
  pal <- c('#ce8600','#3c65eb','#e057ac','#e81988','tan3','yellow','lightblue1',
                                                   '#59be5a','#34b14e','#daf399','#1f8342','#874300','#fdac34','#a79d73','red3','#d08ece')
  for (i in 1:length(nodes2)) {
    node_cols[[nodes2[i]]] <- list(col = pal[i], textcex = label.cex,
                                   textcol = label.col, srt = lab.srt)
  }
  # create nodes and riverplot object
  nodes <- list(nodes1, nodes_middle, nodes2)
  node.limit <- max(unlist(lapply(nodes, length)))
  
  node_Ys <- lapply(1:length(nodes), function(i) {
    seq(1, node.limit, by = node.limit / length(nodes[[i]]))
  })

rp <- makeRiver(c(nodes1, nodes_middle, nodes2), edge_list,
                  node_xpos = node_Xs, node_ypos = unlist(node_Ys), node_styles = node_cols
  )

river.yscale=4

river.lty =0.1

river.node_margin=0.01
p<- riverplot(rp,river.lty=river.lty,yscale=river.yscale,node_margin=river.node_margin)
print(p)

