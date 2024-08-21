library(Seurat)

library(ggplot2)


Idents(drg)<-'CellType'

sub <- subset(drg,idents = c('sensory neuron progenitor',
                             'glia cell progenitor','precursor1','precursor2',
                             'neural crest cell'))

#import auc matrix from regulon analysis

auc_mtx <- read.csv(file='DRG/scenic/human_major_type_relax_auc_mtx.csv',sep='\t')

dim(auc_mtx)

auc_mtx <- t(auc_mtx)

rownames(auc_mtx)<- gsub('\\...','',rownames(auc_mtx))


sub[['regulon']] <- CreateAssayObject(count = as.matrix(auc_mtx[,colnames(sub)]))


tf_target = read.csv(file='DRG/scenic/human_major_type_tf_target.csv',sep='\t')

library(dplyr)
human_deg <- read.csv(file='DRG/HumanMouse/human_drg_major_type_deg_major.csv')
human_deg <- arrange(human_deg,avg_log2FC,.by_group = T)
human_deg <- human_deg%>%group_by(cluster)%>%arrange(desc(avg_log2FC),.by_group = T)

top50_h <-human_deg%>%group_by(cluster)%>%top_n(250,avg_log2FC)

pos <- grep('sensory neuron progenitor|glia cell progenitor|precursor1|precursor2|neural crest',top50_h$cluster)

top50_h <- top50_h[pos,]

Idents(drg1)<-'CellType3'

sub <- subset(drg1,idents = c('sensory neuron progenitor',
                             'glia cell progenitor','precursor1','precursor2',
                             'neural crest cell'))

mat_re = sub@assays$regulon@counts

mat_re <- t(mat_re)

#constructed correlation matrix
re_cor = sparse_cor(mat_re)

pca_mat_re <- irlba::prcomp_irlba(re_cor, n=50)$x
        rownames(pca_mat_re) <- rownames(re_cor)
        x <- as.matrix(pca_mat_re)

#kmeans analysis
cl = kmeans(pca_mat_re, 10, nstart = 5)

g_umap$kmean = as.character(cl$cluster)

g_umap <- umap_tbl <- uwot::umap(pca_mat_re[,c(1:25)]) %>%
        {colnames(.) <- c('UMAP_1', 'UMAP_2'); .} %>%
        as_tibble(rownames='gene')
    


ggplot(g_umap,aes(x=UMAP_1,y=UMAP_2,color=kmean))+geom_point()+theme_bw()+geom_text_repel(aes(label=kmean))

sub <- NormalizeData(sub)

For each k-mean clusters, perform gene-set enrichment analysis
for(i in 1:length(unique(g_umap$kmean))){

sub <- AddModuleScore(sub,features = list(g_umap[g_umap$kmean==i,]$gene),name = paste0('group','i'))

}

#visualization of enrichment  with featureplot 
FeaturePlot(sub,features = 'group1')
FeaturePlot(sub,features = 'group2')
FeaturePlot(sub,features = 'group3')
FeaturePlot(sub,features = 'group4')
FeaturePlot(sub,features = 'group5')


library(reshape2)

g_cor_sub_l = melt(re_cor)


g_cor_sub_l <- data.frame(tf = g_cor_sub_l$variable,target=g_cor_sub_l$target,corr=g_cor_sub_l$value)

#define edge to show
pos <- grep('TRUE',abs(g_cor_sub_l$corr)>0.15)

g_cor_sub_l$edge_type = 'no'
g_cor_sub_l$edge_type[pos] = 'yes'

table(g_cor_sub_l$edge_type)

#define gene to label

pos <- which(g_cor_sub_l$target%in%gene_selected$V1)
length(pos)

g_cor_sub_l$gene_label[pos] = g_cor_sub_l$target[pos]


pos <- which(g_umap$gene%in%gene_selected$V1)
length(pos)

g_umap$gene_label = ''
g_umap$gene_label[pos] = g_umap$gene[pos]


#compute gene average expression
sub$gene_avg = 'drg'

avg_expr <- AverageExpression(sub,group.by = 'gene_avg',features = unique(g_umap$gene),assays = 'RNA')

sub$time = sub$GW

sub$time <- gsub('GW','',sub$time)

sub$time = as.numeric(sub$time)

avg_expr_time <- AverageExpression(sub,group.by = 'time',features = unique(g_umap$gene),assays = 'RNA')


#network construction and visualization
gene_graph <- as_tbl_graph(g_cor_sub_l)%>%activate(edges)%>%mutate(from_node=.N()$name[from], to_node=.N()$name[to],edge_type=edge_type)%>%
activate(nodes)%>%mutate(centrality = centrality_pagerank(),betweeness = centrality_betweenness(),degree=centrality_degree())%>%inner_join(g_umap,by = c('name'='gene'))


p <- ggraph(gene_graph, x=UMAP_1, y=UMAP_2)


p <- p + geom_edge_diagonal(aes(color=edge_type_p1),width=0.25,lineend = 'round') +
            scale_edge_color_manual(values=c('yes'='#DCE5E5','no'='transparent','yes_p1'='grey50'))+
scale_edge_alpha_manual(values=c('yes'=0.1,'no'=0.1,'yes_p1'='1'))

p <- p + geom_edge_diagonal(aes(color=edge_type_p1),width=0.25,lineend = 'round') +
            scale_edge_color_manual(values=c('yes'='#DCE5E5','no'='transparent','yes_p1'='grey50'))+
scale_edge_alpha_manual(values=c('yes'=0.1,'no'=0.1,'yes_p1'='1'))


p <- p + geom_edge_diagonal(aes(color=edge_type),width=0.25,lineend = 'round') +
            scale_edge_color_manual(values=c('yes'='grey70','no'='transparent'))


node_color = pals::magma(100)

p <- p + geom_node_point(aes(fill=pseudotime, size=node_size), shape=21) +
            scale_fill_gradientn(colors=node_color,limits = c(5,21),oob = scales::squish)


table(g_cor_sub_l$edge_type)

#precursor2
p <- p + geom_node_point(aes(size=node_size,alpha=p1,fill=p1), shape=21) +
scale_fill_manual(values=c('precursor1'='#EA9292','other'='#EA9292','precursor1_tf'='#D22928'))+
scale_alpha_manual(values=c('precursor1'=0.1,'other'=0.1,'precursor1_tf'=1))


#precursor1
p <- p + geom_node_point(aes(size=node_size,alpha=p1,fill=p1), shape=21) +
scale_fill_manual(values=c('precursor1'='#A3D8BF','other'='#A3D8BF','precursor1_tf'='#279D68'))+
scale_alpha_manual(values=c('precursor1'=0.1,'other'=0.1,'precursor1_tf'=1))


#snp
p <- p + geom_node_point(aes(size=node_size,alpha=p1,fill=p1), shape=21) +
scale_fill_manual(values=c('precursor1'='#EDB388','other'='#EDB388','precursor1_tf'='#E86C10'))+
scale_alpha_manual(values=c('precursor1'=0.1,'other'=0.1,'precursor1_tf'=1))


#glia progenitor
p <- p + geom_node_point(aes(size=node_size,alpha=p1,fill=p1), shape=21) +
scale_fill_manual(values=c('precursor1'='#9CE8F2','other'='#9CE8F2','precursor1_tf'='blue'))+
scale_alpha_manual(values=c('precursor1'=0.1,'other'=0.1,'precursor1_tf'=1))


graph = 'module_graph'
    layout = 'umap'
    edge_width = 0.2
    edge_color = c('-1'='darkgrey', '1'='orange')
    node_color = pals::magma(100)
    node_size = c(1,5)
    text_size = 10
    color_nodes = TRUE
    label_nodes = TRUE
    color_edges = TRUE

#+scale_edge_width(range = c(0, 0.1))
#print(p)
#pdf('/gpfs2/wulab10/DRG/progenitor_gene_module/progenitor_gene_module.pdf',height = 13,width = 12)
p <- p + geom_node_text(
            aes(label=gene_label_p1),
            repel=T, size=text_size/ggplot2::.pt, max.overlaps=99999)
p <- p + scale_size_continuous(range=c(2,25)) +
       theme_void()
print(p)
#dev.off()



table(g_umap$kmean)

ggplot(g_umap,aes(x=UMAP_1,y=UMAP_2,color=kmean))+geom_point()

p <- ggraph(gene_graph, x=UMAP_1, y=UMAP_2)

p <- p + geom_edge_diagonal(aes(color=edge_type),width=0.25,lineend = 'round') +
            scale_edge_color_manual(values=c('yes'='grey70','no'='transparent'))

p <- p + geom_node_point(aes(fill=kmean, size=node_size), shape=21) +
            scale_fill_manual(values=c('1'='#FF7491','2'='#C3E2AA','3'='#C563FF','4'='#77A799','5'='#FEFF30',
                                       '6'='#41A3A3','7'='#D8165C','8'='#3DBCE0','9'='#D1AF2E',
                                       '10'='#0B74B1'))


p <- p + scale_size_continuous(range=c(2,8)) +
       theme_void()




sub<- AddModuleScore(sub,features = list(g_umap[g_umap$kmean==1,]$gene),name='G1')

sub<- AddModuleScore(sub,features = list(g_umap[g_umap$kmean==2,]$gene),name='G2')

sub<- AddModuleScore(sub,features = list(g_umap[g_umap$kmean==3,]$gene),name='G3')

sub<- AddModuleScore(sub,features = list(g_umap[g_umap$kmean==4,]$gene),name='G4')

sub<- AddModuleScore(sub,features = list(g_umap[g_umap$kmean==5,]$gene),name='G5')

sub<- AddModuleScore(sub,features = list(g_umap[g_umap$kmean==6,]$gene),name='G6')

sub<- AddModuleScore(sub,features = list(g_umap[g_umap$kmean==7,]$gene),name='G7')

sub<- AddModuleScore(sub,features = list(g_umap[g_umap$kmean==8,]$gene),name='G8')

sub<- AddModuleScore(sub,features = list(g_umap[g_umap$kmean==9,]$gene),name='G9')

sub<- AddModuleScore(sub,features = list(g_umap[g_umap$kmean==10,]$gene),name='G10')




levels(sub)<- as.factor(rev(c('neural crest cell','glia cell progenitor','sensory neuron progenitor','precursor1','precursor2')))



p_wide = p$data

p_wide = p_wide[,c(1,3,4)]

colnames(p_wide)<- c('avg.expr','features.plot','subtype')

p_wide$features.plot = as.character(p_wide$features.plot)
p_wide$subtype = as.character(p_wide$subtype)

p_wide = reshape(p_wide, idvar = "features.plot", timevar = "subtype",direction = 'wide')

p_wide <- p_wide[,-1]

library(pheatmap)

100,units='mm',compression = 'lzw',res=600)
pheatmap(p_wide,cluster_cols = F,cluster_rows = F,scale='row',color = colorRampPalette(c('white','grey','black'))(10))

p_wide = melt(p_wide)

p <- p + geom_node_text(
            aes(label=gene_label),
            repel=T, size=text_size/ggplot2::.pt, max.overlaps=99999)
p <- p + scale_size_continuous(range=c(2,8)) +
       theme_void()
print(p)

p <- p + geom_node_text(
            aes(label=gene_label),
            repel=T, size=text_size/ggplot2::.pt, max.overlaps=99999)
p <- p + scale_size_continuous(range=c(2,8)) +
       theme_void()
print(p)

#+scale_edge_width(range = c(0, 0.1))
#print(p)
p <- p + geom_node_text(
            aes(label=gene_label),
            repel=T, size=text_size/ggplot2::.pt, max.overlaps=99999)
p <- p + scale_size_continuous(range=c(2,8)) +
       theme_void()
print(p)

#+scale_edge_width(range = c(0, 0.1))
#print(p)
p <- p + geom_node_text(
            aes(label=gene_label),
            repel=T, size=text_size/ggplot2::.pt, max.overlaps=99999)
p <- p + scale_size_continuous(range=c(2,8)) +
       theme_void()
print(p)



g_umap$p1 = 'other'
pos <- which(g_umap$gene%in%type_marker[type_marker$cluster=='sensory neuron progenitor',]$gene)

g_umap$p1[pos]<-'precursor1'

pos <- which(g_umap[g_umap$p1=='precursor1',]$gene%in%tf_target$tf)
gene = g_umap[g_umap$p1=='precursor1',]$gene[pos]

g_umap[g_umap$p1=='precursor1',]$p1[pos] = 'precursor1_tf'

table(g_umap$p1)

pos <- which(g_cor_sub_l[g_cor_sub_l$tf%in%gene,]$target%in%type_marker[type_marker$cluster=='sensory neuron progenitor',]$gene)


g_cor_sub_l$edge_type_p1 = g_cor_sub_l$edge_type


g_cor_sub_l[g_cor_sub_l$tf%in%gene,]$edge_type_p1 <-'yes_p1'

pos <- grep('no',g_cor_sub_l$edge_type)
g_cor_sub_l$edge_type_p1[pos]<-'no'

g_umap$gene_label_p1 = ''

g_umap[g_umap$p1=='precursor1_tf',]$gene_label_p1 = g_umap[g_umap$p1=='precursor1_tf',]$gene


table(g_umap$gene_label_p1)

Idents(sub)<-'CellType3'

type_marker<- FindAllMarkers(sub,only.pos = T)

pos <-which(g_umap$gene%in%type_marker[type_marker$cluster=='precursor1',]$gene)
length(pos)
g_umap_p1 = g_umap[pos,]

pos <-which(g_cor_sub_l$target%in%type_marker[type_marker$cluster=='precursor1',]$gene)
length(pos)
g_cor_sub_l_p1 =g_cor_sub_l[pos,]

pos <- grep('TRUE',g_cor_sub_l_p1$corr>=0.15)
g_cor_sub_l_p1=  g_cor_sub_l_p1[pos,]

pos <- which(g_umap_p1$gene%in%g_cor_sub_l_p1$target)
length(pos)
g_umap_p1 <- g_umap_p1[pos,]

gene_graph <- as_tbl_graph(g_cor_sub_l_p1)%>%activate(edges)%>%mutate(from_node=.N()$name[from], to_node=.N()$name[to],edge_type=edge_type)%>%
activate(nodes)%>%mutate(centrality = centrality_pagerank(),betweeness = centrality_betweenness(),degree=centrality_degree())%>%inner_join(g_umap_p1,by = c('name'='gene'))


p <- ggraph(gene_graph, x=UMAP_1, y=UMAP_2)


p <- p + geom_edge_diagonal(aes(color=edge_type),width=0.15,lineend = 'round') +
            scale_edge_color_manual(values=c('yes'='grey70','no'='transparent'))


node_color = pals::magma(100)

p <- p + geom_node_point(aes(fill=pseudotime, size=node_size), shape=21) +
            scale_fill_gradientn(colors=node_color,limits = c(5,21),oob = scales::squish)


graph = 'module_graph'
    layout = 'umap'
    edge_width = 0.2
    edge_color = c('-1'='darkgrey', '1'='orange')
    node_color = pals::magma(100)
    node_size = c(1,5)
    text_size = 10
    color_nodes = TRUE
    label_nodes = TRUE
    color_edges = TRUE

#+scale_edge_width(range = c(0, 0.1))
#print(p)
p <- p + geom_node_text(
            aes(label=gene_label),
            repel=T, size=text_size/ggplot2::.pt, max.overlaps=99999)
p <- p + scale_size_continuous(range=c(2,8)) +
       theme_void()
print(p)

pos <-which(g_umap$gene%in%type_marker[type_marker$cluster=='precursor2',]$gene)
length(pos)
g_umap_p2 = g_umap[pos,]

pos <-which(g_cor_sub_l$target%in%type_marker[type_marker$cluster=='precursor2',]$gene)
length(pos)
g_cor_sub_l_p2 =g_cor_sub_l[pos,]

pos <- grep('TRUE',g_cor_sub_l_p2$corr>=0.15)
g_cor_sub_l_p2=  g_cor_sub_l_p2[pos,]


pos <- which(g_umap_p2$gene%in%g_cor_sub_l_p2$target)
length(pos)
g_umap_p2<- g_umap_p2[pos,]

gene_graph <- as_tbl_graph(g_cor_sub_l_p2)%>%activate(edges)%>%mutate(from_node=.N()$name[from], to_node=.N()$name[to],edge_type=edge_type)%>%
activate(nodes)%>%mutate(centrality = centrality_pagerank(),betweeness = centrality_betweenness(),degree=centrality_degree())%>%inner_join(g_umap_p2,by = c('name'='gene'))


p <- ggraph(gene_graph, x=UMAP_1, y=UMAP_2)


p <- p + geom_edge_diagonal(aes(color=edge_type),width=0.15,lineend = 'round') +
            scale_edge_color_manual(values=c('yes'='grey70','no'='transparent'))


node_color = pals::magma(100)

p <- p + geom_node_point(aes(fill=pseudotime, size=node_size), shape=21) +
            scale_fill_gradientn(colors=node_color,limits = c(5,21),oob = scales::squish)


graph = 'module_graph'
    layout = 'umap'
    edge_width = 0.2
    edge_color = c('-1'='darkgrey', '1'='orange')
    node_color = pals::magma(100)
    node_size = c(1,5)
    text_size = 10
    color_nodes = TRUE
    label_nodes = TRUE
    color_edges = TRUE

#+scale_edge_width(range = c(0, 0.1))
#print(p)
p <- p + geom_node_text(
            aes(label=gene_label),
            repel=T, size=text_size/ggplot2::.pt, max.overlaps=99999)
p <- p + scale_size_continuous(range=c(2,8)) +
       theme_void()
print(p)

library(VGAM)

sub$time= sub$GW
sub$time <- gsub('GW','',sub$time)
sub$time = as.numeric(sub$time)

responseMatrix <- function(models, newdata = NULL, response_type="response", cores = 1) {
  res_list <- mclapply(models, function(x) {
    if (is.null(x)) { NA } else {
      if (x@family@vfamily %in% c("negbinomial", "negbinomial.size")) {
        predict(x, newdata = newdata, type = response_type)
      } else if (x@family@vfamily %in% c("uninormal")) {
        predict(x, newdata = newdata, type = response_type)
      }
      else {
        10^predict(x, newdata = newdata, type = response_type)
      }
    }
  }, mc.cores = cores)
  
  res_list_lengths <- lapply(res_list[is.na(res_list) == FALSE],
                             length)
  stopifnot(length(unique(res_list_lengths)) == 1)
  num_na_fits <- length(res_list[is.na(res_list)])
  if (num_na_fits > 0) {
    na_matrix <- matrix(rep(rep(NA, res_list_lengths[[1]]),
                            num_na_fits), nrow = num_na_fits)
    row.names(na_matrix) <- names(res_list[is.na(res_list)])
    non_na_matrix <- Matrix::t(do.call(cbind, lapply(res_list[is.na(res_list) ==
                                                                FALSE], unlist)))
    row.names(non_na_matrix) <- names(res_list[is.na(res_list) ==
                                                 FALSE])
    res_matrix <- rbind(non_na_matrix, na_matrix)
    res_matrix <- res_matrix[names(res_list), ]
  }
  else {
    res_matrix <- Matrix::t(do.call(cbind, lapply(res_list, unlist)))
    row.names(res_matrix) <- names(res_list[is.na(res_list) ==
                                              FALSE])
  }
  res_matrix
}



