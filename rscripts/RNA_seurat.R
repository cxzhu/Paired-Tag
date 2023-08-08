
###
library(Seurat)
setwd("/projects/ps-renlab/chz272/transposon/05.Paired-ChIP/20.NovaSeq/04.Single_cell_All/04.RNA_Proc/01.RNA_seurat")
#pt.seu<-Read10X(data.dir="/projects/ren-transposon/home/chz272/transposon/05.Paired-ChIP/20.NovaSeq/04.Single_cell_All/03.filtered_matrices/RNA_filtered_matrix")
#pt.seu<-Read10X(data.dir="/projects/ren-transposon/home/chz272/transposon/05.Paired-ChIP/20.NovaSeq/04.Single_cell_All/03.filtered_matrices/RNA_filtered_matrix")
#pt.s<-CreateSeuratObject(counts=pt.seu, project="Paired-tag RNA", min.cells=0, min.features=0)
pt.s

pt.s[["percent.mt"]] <- PercentageFeatureSet(pt.s, pattern = "^mt-")
pt.s<-NormalizeData(pt.s, normalization.method="LogNormalize", scale.factor=10000)
pt.s<-FindVariableFeatures(pt.s, selection.method="vst", nfeatures=25000)
pt.s <- ScaleData(pt.s, vars.to.regress = "percent.mt")


pt.s<-FindVariableFeatures(pt.s, selection.method="vst", nfeatures=25000)
pt.s <- RunPCA(pt.s, features = VariableFeatures(object = pt.s))
ElbowPlot(pt.s, ndims = 50)
pc.use<-c(1:27)
k=15
pt.s <- FindNeighbors(pt.s, dims = pc.use, k.param=k)
pt.s <- FindClusters(pt.s, resolution = 1.5, algorithm = 1)
#pt.s <- RunUMAP(pt.s, dims = pc.use, metric = "euclidean", min.dist = 0.01, local.connectivity=10, spread=1, seed.use=131, umap.method = "umap-learn", n.neighbors = k, uwot.sgd=TRUE, verbose=TRUE)
DimPlot(pt.s, reduction = "umap", label=T, pt.size=0.1)
plot(pt.s@reductions$umap@cell.embeddings, pch=19, cex=0.01, col=col.batch[id.batch])


output<-(pt.s@reductions$umap@cell.embeddings)
output<-cbind(output, pt.s@meta.data$seurat_clusters)
colnames(output)<-c("UMAP_1","UMAP_2", "Membership")
output[,3]<-output[,3]-1
write.table(output, file="2020_03_31_RNA_UMAP_37.xls", sep="\t", quote=F, row.names=T, col.names=T)
#saveRDS(pt.s, file="2020_03_31_RNA_Seurat.rds")



rna.meta<-read.csv("/projects/ren-transposon/home/chz272/transposon/05.Paired-ChIP/20.NovaSeq/04.Single_cell_All/04.RNA_Proc/01.RNA_seurat/2020_04_01_RNA_summary.xls", sep="\t", head=T)
batch<-rep(1,dim(rna.meta)[1])
batch[rna.meta$Target=="H3K27me3"]<-2
batch[rna.meta$Target=="H3K9me3"]<-2
pt.s@meta.data$batch<-batch
library(harmony)
options(repr.plot.height = 5, repr.plot.width = 12)
p1<-DimPlot(object=pt.s, reduction="pca", pt.size=.1, group.by="batch") #, do.return = TRUE)
p2<-VlnPlot(object=pt.s, features="PC_1", group.by="batch", pt.size=.1)
plot_grid(p1,p2)

options(repr.plot.height = 2.5, repr.plot.width = 6)
pt.s<-pt.s %>%
  RunHarmony("batch", plot_convergence=TRUE)
harmony_embeddings<-Embeddings(pt.s, 'harmony')
harmony_embeddings[1:5, 1:5]

ElbowPlot(pt.s,reduction="harmony", ndims=50)

pc.use<-c(1:27)
pt.s<-pt.s %>%
  RunUMAP(reduction = "harmony", dims=pc.use, metric = "euclidean", min.dist = 0.01, n.neighbors=k, uwot.sgd=TRUE, verbose=TRUE, local.connectivity=10, spread=1, seed.use=131) %>%
  FindNeighbors(reduction = "harmony", dims=pc.use, k.param=k) %>%
  FindClusters(resolution = 1, algorithm=2 ) %>%
  identity()


k=15
pc.use<-1:27

pt.s <- FindNeighbors(pt.s,reduction="harmony", dims=pc.use, k.param=k)
pt.s <- FindClusters(pt.s, resolution = 1.5, algorithm = 2)

pt.s <- RunUMAP(pt.s, dims = pc.use, metric = "euclidean", min.dist = 0.01, local.connectivity=10, spread=1, seed.use=131, umap.method = "uwot", n.neighbors = k, uwot.sgd=TRUE, verbose=TRUE)


options(repr.plot.height=4, repr.plot.width=10);DimPlot(pt.s, reduction = "umap", group.by= "batch", pt.size= .1, split.by="batch")
options(repr.plot.height=4, repr.plot.width=10);DimPlot(pt.s, reduction = "umap", group.by= "origin", pt.size= .1, split.by="origin")
options(repr.plot.height=4, repr.plot.width=6);DimPlot(pt.s, reduction="umap", label=TRUE, pt.size=.1)

