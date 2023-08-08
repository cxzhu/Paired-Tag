library(RColorBrewer)
library(Matrix)
col2<-brewer.pal(n=12, name="Paired")
col1<-brewer.pal(n=9, name="Set1")
col<-c(col2,col1)

plot_cell_ID<-function(data, cutoff=100, sample="NA", ...){
  plot(rev(sort(data)), pch=19, cex=0.5, col="grey", xlab="# of Cellular Barcodes", ylab="# of Reads", log="xy", main=sample)
  filt<-data[data>cutoff];med=median(data[data>cutoff])
  points(rev(sort(filt)), pch=19, cex=0.5, col="red")
  lines(c(1, length(data)), c(med,med), lty=2, col="grey", lwd=1)
  frac<-sum(data[data>cutoff])/sum(data); frac<-(as.integer(frac*10000)); frac<-frac/100
  legend("bottomleft", legend=c( paste("Total:", length(data), sep="\t"), paste("PF:", length(filt), sep="\t\t"), paste("Reads:\t", frac, "%", sep=""), paste("Median:\t", med, sep="")), pch=19, col=c("grey", "red", "white", "white"), bty="n", cex=1)
}

plot_cell_ID_rna<-function(data, cutoff=100){
  data<-data*3
  
  plot(rev(sort(data)), pch=19, cex=0.5, col="grey", xlab="# of Cellular Barcodes", ylab="# of Reads", log="xy", ylim=c(1,100000))
  filt<-data[data>cutoff];med=median(data[data>cutoff])
  points(rev(sort(filt)), pch=19, cex=0.5, col="red")
  lines(c(1, length(data)), c(med,med), lty=2, col="grey", lwd=1)
  frac<-sum(data[data>cutoff])/sum(data); frac<-(as.integer(frac*10000)); frac<-frac/100
  legend("bottomleft", legend=c( paste("Total:", length(data), sep="\t"), paste("PF:", as.integer(length(filt)), sep="\t\t"), paste("Reads:\t", frac, "%", sep=""), paste("Median:\t", med, sep="")), pch=19, col=c("grey", "red", "white", "white"), bty="n")
}

plot_BC_col<-function(data, cutoff, sample){
  data<-data[(data[,2]+data[,3])>cutoff,]
  xlab="UMI hg19";ylab="UMI mm10"
  xlim<-c(0, max(c(data[,2], data[,3])));ylim<-xlim
  human<- data[,2] > data[,3]*4
  mouse<- data[,2]*4 < data[,3]
  mix <- (data[,2]*4 > data[,3]) & (data[,2] < data[,3]*4)
  col<-rep("grey", dim(data)[1])
  col[human]<-"blue";col[mouse]<-"red"
  nHuman<-length(human[human]);nMouse<-length(mouse[mouse]);nMix<-length(mix[mix])
  max<-max(max(data[,2]), max(data[,3]))
  plot(data[,2], data[,3], pch=19, cex=0.5, xlim=c(0,max/1), ylim=c(0,max/1), col=col, main=sample, xlab="# of hg19 mapped reads", ylab="# of mm10 mapped reads")
  nMix<-as.integer(nMix/1)
  legend("topright", legend=c( paste("Human cell:", nHuman, sep=" "), paste("Mouse cell:", nMouse, sep=" "), paste("Mixed cell:", nMix, sep= " ") ), pch=19, col=c("blue", "red", "grey"))
}


plot_BC_col_RNA<-function(data, cutoff.human, cutoff.mouse,cutoff2, sample){
  data<-data[data[,2]>cutoff.human | data[,3]>cutoff.mouse,]
  data<-data[(data[,2]+data[,3])<cutoff2,]
  xlab="UMI hg19";ylab="UMI mm10"
  xlim<-c(0, max(c(data[,2], data[,3])));ylim<-xlim
  human<- data[,2] > data[,3]*4 & data[,2]>cutoff.human
  mouse<- data[,2]*4 < data[,3] & data[,3]>cutoff.mouse
  mix <- (data[,2]*4 > data[,3]) & (data[,2] < data[,3]*4) &  data[,2]>cutoff.human & data[,3]>cutoff.mouse
  col<-rep("white", dim(data)[1])
  col[human]<-"blue";col[mouse]<-"red";col[mix]<-"grey"
  nHuman<-length(human[human]);nMouse<-length(mouse[mouse]);nMix<-length(mix[mix])
  max<-max(c(data[,2], data[,3]))*1.25
  plot(data[,2], data[,3], pch=19, cex=0.5, xlim=c(0,max), ylim=c(0,max), col=col, main=sample, xlab="# of hg19 mapped reads", ylab="# of mm10 mapped reads")
  nMix<-as.integer(nMix/1)
  legend("topright", legend=c( paste("Human cell:", nHuman, sep=" "), paste("Mouse cell:", nMouse, sep=" "), paste("Mixed cell:", nMix, sep= " ") ), bty="n",pch=19, col=c("blue", "red", "grey"))
}

plot_BC_col_K4me1<-function(data, cutoff, sample){
  data<-data[(data[,2]+data[,3])>cutoff,]
  xlab="UMI hg19";ylab="UMI mm10"
  xlim<-c(0, max(c(data[,2], data[,3])));ylim<-xlim
  human<- data[,2] > data[,3]*4
  mouse<- data[,2]*4 < data[,3]
  mix <- (data[,2]*4 > data[,3]) & (data[,2] < data[,3]*4)
  col<-rep("grey", dim(data)[1])
  col[human]<-"blue";col[mouse]<-"red"
  nHuman<-length(human[human]);nMouse<-length(mouse[mouse]);nMix<-length(mix[mix])
  max<-max(max(data[,2]), max(data[,3]))
  plot(data[,2], data[,3], pch=19, cex=0.5, xlim=c(0,max/2), ylim=c(0,max/2), col=col, main=sample, xlab="# of hg19 mapped reads", ylab="# of mm10 mapped reads")
  nMix<-as.integer(nMix/1)
  legend("topright", legend=c( paste("Human cell:", nHuman*2, sep=" "), paste("Mouse cell:", nMouse*2, sep=" "), paste("Mixed cell:", nMix, sep= " ") ), pch=19, col=c("blue", "red", "grey"))
}
plot_BC_col_K27ac<-function(data, cutoff, sample){
  data<-data[(data[,2]+data[,3])>cutoff,]
  xlab="UMI hg19";ylab="UMI mm10"
  xlim<-c(0, max(c(data[,2], data[,3])));ylim<-xlim
  human<- data[,2] > data[,3]*4
  mouse<- data[,2]*4 < data[,3]
  mix <- (data[,2]*4 > data[,3]) & (data[,2] < data[,3]*4)
  col<-rep("grey", dim(data)[1])
  col[human]<-"blue";col[mouse]<-"red"
  nHuman<-length(human[human]);nMouse<-length(mouse[mouse]);nMix<-length(mix[mix])
  max<-max(max(data[,2]), max(data[,3]))
  plot(data[,2], data[,3], pch=19, cex=0.5, xlim=c(0,max/1), ylim=c(0,max/1), col=col, main=sample, xlab="# of hg19 mapped reads", ylab="# of mm10 mapped reads")
  nMix<-as.integer(nMix/1)
  legend("topright", legend=c( paste("Human cell:", nHuman*3, sep=" "), paste("Mouse cell:", nMouse*2, sep=" "), paste("Mixed cell:", nMix, sep= " ") ), pch=19, col=c("blue", "red", "grey"))
}

test<-function(){
  plot(c(0,0), xlim=c(0,100), ylim=c(0,100), cex=0.5, pch=19, col="grey")
  lines(c(0,100), c(10,10), lty=2)
  
}


plot_BC<-function(dna, rna, c1=300, c2=500, title=title){
  layout(matrix(c(1,2,1,3),2,2,byrow=TRUE))
  m<-merge(dna, rna, by="V1")
  plot(m[,2:3],pch=19, cex=0.1, log="xy", xlab="DNA reads", ylab="RNA reads", col="grey", xlim=c(1,100000), ylim=c(1,100000), main=title)
  lines(c(c1,c1),c(1,5000000), lty=2);lines(c(1,5000000),c(c2,c2), lty=2)
  pf<-m[m[,2] >= c1 & m[,3] >= c2,]
  points(pf[,2:3], pch=19, cex=0.1, col="firebrick")
  legend("bottomright", legend=c( paste("DNA cutoff:", c1, sep=" "), paste("RNA cutoff:", c2, sep = " "), paste("# PF:", dim(pf)[1], sep=" ")), bty="n")
  data<-na.omit(pf)
  dna<-log10(data[,2])
  rna<-log10(data[,3])
  hist(dna, xlim=c(1,6), col="deepskyblue", breaks=50, border=F, main="DNA # of loci per nucle", xlab="# of loci, log10")
  m.dna<-median(dna)
  m.dna.v<-median(data[,2]);n.cell=dim(data)[1]
  lines(c(m.dna,m.dna),c(0,5000), lty=2, lwd=2, col="black")
  legend("topright", legend=c(paste("Median #:", m.dna.v, sep=" "), paste("# of nuclei:", n.cell, sep=" ")), bty="n")
  m.rna<-median(rna)
  m.rna.v<-median(data[,3]);n.cell=dim(data)[1]
  hist(rna, xlim=c(1,6), col="forestgreen", breaks=50, border=F, main="RNA # of UMI per nucle", xlab="# of UMI, log10")
  lines(c(m.rna,m.rna),c(0,5000), lty=2, lwd=2, col="black")
  legend("topright", legend=c(paste("Median #:", m.rna.v, sep=" "), paste("# of nuclei:", n.cell, sep=" ")), bty="n")
  val.cell<-data[,1]
  write.table(val.cell, file=paste(title,"filt.xls", sep="."), col.names=F, row.names=F, sep="\t", quote=F)
}


plot_individual<-function(r, cutoff.dna=200, cutoff.rna=200, name=""){
  layout(matrix(c(1,1,2,2,1,1,2,2,3,3,4,4,3,3,5,5),4,4,byrow=TRUE))
  f.d<-r[r[,2]>cutoff.dna,];med.d<-median(f.d[,2]);num.d<-dim(f.d)[1];fr.d<-sum(f.d[,2])/sum(r[,2]);fr.d<-as.integer(fr.d*10000);fr.d<-fr.d/100
  f.r<-r[r[,3]>cutoff.rna,];med.r<-median(f.r[,3]);num.r<-dim(f.r)[1];fr.r<-sum(f.d[,3])/sum(r[,3]);fr.r<-as.integer(fr.r*10000);fr.r<-fr.r/100
  f.b<-r[r[,2]>cutoff.dna & r[,3]>cutoff.rna,];med.jd<-median(f.b[,2]);med.jr<-median(f.b[,3]);num.j<-dim(f.b)[1]
  plot(rev(sort(r[,2])), pch=19, cex=0.25, col="grey", log="xy", xlab="Index", ylab="# of DNA UMI");points(rev(sort(f.d[,2])), pch=19, cex=0.25, col="red")
  legend("bottomleft", pch=19, col=c("grey", "red", "white", "white"),legend=c(paste("Total:", dim(r)[1], sep="\t"), paste("PF:", num.d, sep="\t\t"), paste("%Reads:", fr.d, sep="\t"), paste("Median:", med.d, sep="\t")) , bty="n")
  plot(rev(sort(r[,3])), pch=19, cex=0.25, col="grey", log="xy", xlab="Index", ylab="# of RNA UMI");points(rev(sort(f.r[,3])), pch=19, cex=0.25, col="red")
  legend("bottomleft", pch=19, col=c("grey", "red", "white", "white"),legend=c(paste("Total:", dim(r)[1], sep="\t"), paste("PF:", num.r, sep="\t\t"), paste("%Reads:", fr.r, sep="\t"), paste("Median:", med.r, sep="\t")) , bty="n")
  plot(r[,2:3], pch=19, cex=0.25, col="grey", log="xy",xlab="# of DNA UMI", ylab="# of RNA UMI", main=paste(name, "", sep="\t"))
  points(f.b[,2:3], pch=19, cex=0.25, col="red");lines(c(1,1000000), c(cutoff.rna, cutoff.rna), lty=2, lwd=1, col="grey25");lines( c(cutoff.dna, cutoff.dna), c(1,1000000),lty=2, lwd=1, col="grey25")
  legend("bottomright", legend=c(paste("Total:", dim(r)[1], sep="\t"), paste("PF:", dim(f.b)[1], sep="\t\t")), pch=19, col=c("grey", "red"), bty="n")
  hist(log10(f.b[,2]), xlim=c(1,6), col="deepskyblue", breaks=50, border=F, main="DNA # of loci per nucle", xlab="# of loci, log10")
  lines(c(log10(med.jd), log10(med.jd)),c(0,5000), lty=2, lwd=2, col="grey25")
  legend("topright", legend=c(paste("Median #:", med.jd, sep=" ")), bty="n")
  hist(log10(f.b[,3]), xlim=c(1,6), col="forestgreen", breaks=50, border=F, main="RNA # of loci per nucle", xlab="# of loci, log10")
  lines(c(log10(med.jr), log10(med.jr)),c(0,5000), lty=2, lwd=2, col="grey25")
  legend("topright", legend=c(paste("Median #:", med.jr, sep=" ")), bty="n")
  colnames(f.b)<-c("Cell_ID", "DNA_Loci", "RNA_UMI")
  write.table(f.b, col.names=F, row.names=F, sep="\t", quote=F, file=paste(name,  "filtBC.xls", sep="_"))
}

qc_matrix<-function(dna, rna, cutoff.dna=200, cutoff.rna=200, name=""){
  dna.m<-readMM(paste(dna,"/matrix.mtx", sep=""))
  rna.m<-readMM(paste(rna, "/matrix.mtx", sep=""))
  dna.b<-read.csv(paste(dna,"/barcodes.tsv",sep=""), sep="\t", head=F)
  rna.b<-read.csv(paste(rna,"/barcodes.tsv",sep=""), sep="\t", head=F)
  dna.s<-colSums(dna.m);rna.s<-colSums(rna.m);names(dna.s)<-dna.b[,1];names(rna.s)<-rna.b[,1]
  merge<-merge(dna.s,rna.s,by=0)
  plot_individual(merge, cutoff.dna, cutoff.rna, name)
  return(merge)
}

write_matrix<-function(dna, rna, cutoff.dna=200, cutoff.rna=200){
  dna.m<-readMM(paste(dna,"/matrix.mtx", sep=""))
  rna.m<-readMM(paste(rna, "/matrix.mtx", sep=""))
  dna.b<-read.csv(paste(dna,"/barcodes.tsv",sep=""), sep="\t", head=F)
  rna.b<-read.csv(paste(rna,"/barcodes.tsv",sep=""), sep="\t", head=F)
  dna.c<-read.csv(paste(dna,"/genes.tsv",sep=""), sep=" ", head=F)
  rna.c<-read.csv(paste(rna,"/genes.tsv",sep=""), sep=" ", head=F)
  dna.s<-colSums(dna.m);rna.s<-colSums(rna.m);names(dna.s)<-dna.b[,1];names(rna.s)<-rna.b[,1]
  rownames(dna.m)<-dna.c[,1];rownames(rna.m)<-rna.c[,1]
  colnames(dna.m)<-dna.b[,1];colnames(rna.m)<-rna.b[,1]
  merge<-merge(dna.s,rna.s,by=0)
  rownames(merge)<-merge[,1]
  merge<-merge[merge[,2]>cutoff.dna & merge[,3]>cutoff.rna,]
  dna.f<-dna.m[,rownames(merge)];rna.f<-rna.m[,rownames(merge)]
  s<-rowSums(dna.f);dna.f<-dna.f[s>1,]
  s<-rowSums(rna.f);rna.f<-rna.f[s>1,]
  system(paste("mkdir ", dna, "_filtered", sep=""))
  system(paste("mkdir ", rna, "_filtered", sep=""))
  writeMM(dna.f,paste(dna,"_filtered/matrix.mtx", sep=""))
  write.table(rownames(dna.f), file=paste(dna,"_filtered/genes.tsv", sep=""), quote=F, col.names=F, row.names=F, sep="\t")
  write.table(rownames(merge), file=paste(dna,"_filtered/barcodes.tsv", sep=""), quote=F, col.names=F, row.names=F, sep="\t")
  writeMM(rna.f,paste(rna,"_filtered/matrix.mtx", sep=""))
  write.table(rownames(rna.f), file=paste(rna,"_filtered/genes.tsv", sep=""), quote=F, col.names=F, row.names=F, sep="\t")
  write.table(rownames(merge), file=paste(rna,"_filtered/barcodes.tsv", sep=""), quote=F, col.names=F, row.names=F, sep="\t")


}







print("Load cxzhu.R sucessful")



