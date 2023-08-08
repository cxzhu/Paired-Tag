overlap_coefficient<-function(mat){
  mat.return<-mat
  for(m in 1:dim(mat)[1]){
    for(n in 1:dim(mat)[2]){
      r<-mat[m,n]/min(sum(mat[m,]),sum(mat[,n]))
      mat.return[m,n]<-r
    }
  }
  return(mat.return)
}




setwd("/projects/ren-transposon/home/chz272/transposon/05.Paired-ChIP/20.NovaSeq/04.Single_cell_All/05.DNA_proc")
source("/projects/ren-transposon/home/chz272/transposon/R_projects/Paired-map/R/Paired-map.R")

rna.meta<-read.csv("../04.RNA_Proc/01.RNA_seurat/2020_04_03_Harmony_RNA.metadata.xls", sep="\t" ,head=T);rownames(rna.meta)<-rna.meta$Cell_ID
rownames(rna.meta)<-rna.meta[,1]

pt.k4me1<-loadPairedDataSet(dna="../03.filtered_matrices/04.C10/H3K4me1_filtered_matrix/",
                      rna="NA",
                      project="Paired-tag",
                      min.umi=0,
                      min.fragments=0)

cell.pf<-read.csv("/projects/ps-renlab/chz272/05.Paired-tag/07.split_bam/04.split_DNA_according_to_RNA_clusters_filter_highpileup/H3K4me1.xls", head=T, sep="\t")
is.pf<-match(as.vector(cell.pf[,1]), as.vector(pt.k4me1@barcode));pt.k4me1@barcode<-pt.k4me1@barcode[is.pf];pt.k4me1@mat.raw@mat.dna<-pt.k4me1@mat.raw@mat.dna[,is.pf]

pt.k4me1<-filtMatrix(pt.k4me1, input="dna", low.threshold = 0.5, high.threshold = 2, method="binary")
pt.k4me1<-runDistance(pt.k4me1, input="dna", cell.downsample=1, feature.downsample = 1, method="jaccard");gc()
pt.k4me1<-runNormDistance(pt.k4me1, input="dna");gc()
pt.k4me1<-runPCA(pt.k4me1, input="dna")
plotPCA(pt.k4me1, input="dna")

pc1.dna=1;pc2.dna=25;pc.use.dna<-c(pc1.dna:pc2.dna);# pc.use.rna<-c(pc1.rna:13,15:pc2.rna)
pt.k4me1<-runUMAP(pt.k4me1, use.dims=pc.use.dna, k=15, min_dist=0.01, input="dna", scale="none", local_connectivity=100)
plot(pt.k4me1@umap@vis.dna, pch=19, cex=0.1, col=col.type[rna.id], bty="l", xlab="UMAP_1", ylab="UMAP_2");rna.id<-rna.meta[pt.k4me1@barcode,17]

#umap<-pt.k4me1@umap@vis.dna

rna.id<-rna.meta[pt.k4me1@barcode,17]
pt.k4me1@cluster@id.rna=rna.id
plotFeature(pt.k4me1, feature="depth", feature.type="dna", embedding.type="dna")

meta.k4me1<-rna.meta[as.vector(pt.k4me1@barcode),]
meta.k4me1<-cbind(meta.k4me1,pt.k4me1@umap@vis.dna)
colnames(meta.k4me1)<-c(colnames(meta.k4me1)[1:18], "DNA_UMAP_1","DNA_UMAP_2")

#write.table(meta.k4me1, col.names=T, row.names=F, quote=F, sep="\t", file="2020_04_07_k4me1_metadata.xls")
pdf("./02.02.DNA_H3K4me1_Umap.pdf", width=8, height=5)
library(RColorBrewer)
col.fc<-colorRampPalette(c( "deepskyblue4", "deepskyblue1","paleturquoise2", "paleturquoise3", "paleturquoise4"))(7)
col.fc<-c("deepskyblue4", "deepskyblue3", "deepskyblue2", "darkslategray3", "darkslategray4", "darkslategrey", "lightskyblue4")
col.hc<-colorRampPalette(c("olivedrab4","olivedrab2", "springgreen3","springgreen4"))(4)
col.in<-colorRampPalette(c("orchid2","orchid4"))(3)
col.no<-colorRampPalette(c("tomato4", "tomato", "orange","orange4"))(8)
col.type<-c(col.fc,col.hc,col.in,col.no,"grey35", "grey65", rgb(0,0,0,0))
cluster.name<-c("FC_ExNeu_L2/3", "FC_ExNeu_L5a_Rorb", "FC_ExNeu_L5a_Deptor", "FC_ExNeu_L5_Parm", "FC_ExNeu_L5b_Fezf2", "FC_ExNeu_L6_Syt6", "FC_ExNeu_Claustrum", 
                "HC_ExNeu_CA1", "HC_ExNeu_Subiculum", "HC_ExNeu_CA2/3", "HC_ExNeu_DG", 
                "InNeu_Vip", "InNeu_Sst", "InNeu_Pvalb", 
                "OPC", "Oligo_Il33", "Oligo_Man1a", "Astro_Myoc", "Astro_Nnat", "Microglia", "Endothelial", "Chroid_Plexus", "Neurogenesis_Sox4", "Doublets")
par(mar=c(2.5, 2.5, 3, 15),xpd=TRUE, col.axis="black", cex.axis=1, lab=c(3,3,5), mgp=c(1.25,0.25,0))
plot(meta.k4me1[,19:20], pch=19, cex=0.2, col=col.type[meta.k4me1$ID_Converted], bty="l")
legend("right", legend=cluster.name[1:22], pch=19, cex=1, bty="n", inset=c(-0.5,0), ncol=1, col=col.type)
dev.off()

#saveRDS(pt.k4me1, file="2020_04_07_K4me1_pm.rds")
#pt.k4me1<-readRDS("2020_04_07_K4me1_pm.rds")
meta.k4me1<-read.csv("2020_04_07_k4me1_metadata.xls", sep="\t", head=T);rownames(meta.k4me1)<-(meta.k4me1$Cell_ID)

pt.k4me1<-runCluster(pt.k4me1, use.dims=pc.use.dna, k=7, input="dna");max(pt.k4me1@cluster@id.dna)
plotFeature(pt.k4me1, feature="cluster", feature.type="dna")

plot(meta.k4me1[,19:20], pch=19, cex=0.2, col=col.type[meta.k4me1$ID_Converted], bty="l")
par(mfrow=c(1,2))
plot(meta.k4me1[,19:20], pch=19, cex=0.2, col=colPanel[pt.k4me1@cluster@id.dna], bty="l");legend("top", ncol=5,legend=c(1:max(pt.k4me1@cluster@id.dna)),pch=19, cex=0.75,col=colPanel, bty="n")
plot(meta.k4me1[,19:20], pch=19, cex=0.2, col=col.type[meta.k4me1$ID_Converted], bty="l");legend("top", ncol=4, legend=c(1:22), pch=19, cex=0.75, col=col.type, bty="n")
par(mfrow=c(1,1))

plot(meta.k4me1[,19:20], pch=19, cex=0.2, col=col.type[pt.k4me1@cluster@id.dna], bty="l");legend("top", ncol=5,legend=c(1:max(pt.k4me1@cluster@id.dna)),pch=19, cex=0.75,col=colPanel, bty="n")
id.conv.table<-c(9,18,8,11,4,18,3,20,3,1,6,13,13,25,16,15,1,12,18,2,21,18,11,10,18,25,22,2,18,2,18,16,5)
dna.id.cov<-id.conv.table[pt.k4me1@cluster@id.dna]
plot(meta.k4me1[,19:20], pch=19, cex=0.2, col=col.type[dna.id.cov], bty="l");legend("top", ncol=5,legend=c(1:max(dna.id.cov)),pch=19, cex=0.75,col=col.type, bty="n")
label.table<-c("FC_ExNeu_L2/3", "FC_ExNeu_L5a_Rorb", "FC_ExNeu_L5a_Deptor", "FC_ExNeu_L5_Parm", "FC_ExNeu_L5b_Fezf2", "FC_ExNeu_L6_Syt6", "FC_ExNeu_Claustrum", 
               "HC_ExNeu_CA1", "HC_ExNeu_Subiculum", "HC_ExNeu_CA2/3", "HC_ExNeu_DG", 
               "InNeu_Vip", "InNeu_Sst", "InNeu_Pvalb", 
               "OPC", "Oligo_Il33", "Oligo_Man1a", "Astro_Myoc", "Astro_Nnat", "Microglia", "Endothelial", "Chroid_Plexus", "Doublets", "Doublets", "Lowquality")
app<-cbind(dna.id.cov, label.table[dna.id.cov]);colnames(app)<-c("DNA_ID","DNA_Annotation")
meta.k4me1<-cbind(meta.k4me1, app)



