source("/projects/ren-transposon/home/chz272/transposon/R_projects/Paired-map/R/cxzhu.R")
source("/projects/ren-transposon/home/chz272/transposon/R_projects/Paired-map/R/Paired-map.R")
library(Seurat)

setwd("/projects/ren-transposon/home/chz272/transposon/05.Paired-ChIP/20.NovaSeq/16.Integration/2020_10")
meta<-read.csv("/projects/ps-renlab/chz272/05.Paired-tag/11.2020_20_02_new_clusters/2020_10_06_Paired-Tag_MetaData.xls", sep="\t", head=T)
rna<-readRDS("/projects/ps-renlab/chz272/05.Paired-tag/11.2020_20_02_new_clusters/2020_10_02_RNA_seurat.rds")
rna.pca<-rna@reductions$pca@cell.embeddings





### H3K4me1
pt<-loadPairedDataSet(dna="/projects/ren-transposon/home/chz272/transposon/05.Paired-ChIP/20.NovaSeq/04.Single_cell_All/03.filtered_matrices/04.C10/H3K4me1_filtered_matrix",
                      rna="/projects/ren-transposon/home/chz272/transposon/05.Paired-ChIP/20.NovaSeq/04.Single_cell_All/03.filtered_matrices/04.C10/H3K4me1_RNA_filtered_matrix",
                      project="Paired-tag",
                      min.umi=0,max.umi=1000000000,
                      min.genes=0,min.cells.dna=0, min.cells.rna=0,
                      min.fragments=0,max.fragment=1000000000)


pt.bk<-pt

meta.k4me1<-meta[meta$Target=="H3K4me1",];meta.k4me1<-na.omit(meta.k4me1)
is.pf<-match(as.vector(meta.k4me1$Cell_ID), as.vector(pt@barcode));pt@barcode<-pt@barcode[is.pf];pt@mat.raw@mat.dna<-pt@mat.raw@mat.dna[,is.pf];pt@mat.raw@mat.rna<-pt@mat.raw@mat.rna[,is.pf]


pt<-filtMatrix(pt, input="rna", low.threshold = -0.5, high.threshold = 2.5, method="logMedian")
pt<-runDistance(pt, input="rna", cell.downsample=1, feature.downsample = 1, method="euclidean");gc()
pt<-runNormDistance(pt, input="rna", method="none");gc()

pt<-filtMatrix(pt, input="dna", low.threshold = 0.5, high.threshold = 2, method="binary")
pt<-runDistance(pt, input="dna", cell.downsample=1, feature.downsample = 1, method="jaccard");gc()
pt<-runNormDistance(pt, input="dna");gc()

pt<-runPCA(pt, input="dna")
plotPCA(pt, input="dna")

pt<-runPCA(pt, input="rna")
plotPCA(pt, input="rna")
pc.use<-c(1:25)
pt<-runCluster(pt, use.dims=pc.use, k=15, input="rna");max(pt@cluster@id.rna)
pt<-runUMAP(pt, use.dims=pc.use, k=15, min_dist=0.1, input="rna", scale="none")

plotFeature(pt, feature="cluster", feature.type="rna", embedding.type="rna")
plot(pt@umap@vis.rna, pch=19, cex=0.2, col=col.type[meta.k4me1$ID_Converted], bty="l")


pt<-runFusionMatrices(pt, use.rna=T, use.dna=T, use.tag=F, method="hadamard")



pt<-runPCA(pt, input="int")
plotPCA(pt, input="int")
pc.use<-c(1:50)
pt<-runCluster(pt, use.dims=pc.use, k=15, input="int");max(pt@cluster@id.int)

library(harmony)


pt<-runUMAP(pt, use.dims=pc.use, k=15, min_dist=0.1, input="int", scale="none")


plotFeature(pt, feature="cluster", feature.type="int", embedding.type="int")
plotFeature(pt, feature="Pdgfra", feature.type="rna", embedding.type="int")
plotFeature(pt, feature="Mbp", feature.type="rna", embedding.type="int")



