#### Paired-map
#### Paired-seq/tag associated Multi-modal Analysis Pipeline
#### runJaccard and normOVE borrowed from SnapATAC [https://github.com/r3fang/SnapATAC]
#### Chenxu Zhu, cxzhu@pku.edu.cn
#### 2020-01-29

version="2020.02.05"

###
{ ### External packages
	print("Checking, installing and loading packages...")
	if (!require(Matrix)) install.packages('Matrix')
	library(Matrix)
	if (!require(irlba)) install.packages('irlba')
	library(irlba)
	if (!require(uwot)) install.packages('uwot')
	library(uwot)
	if (!require(RColorBrewer)) install.packages('RColorBrewer')
	library(RColorBrewer)
	if (!require(Rtsne)) install.packages('Rtsne')
	library(Rtsne)
	if (!require(stats)) install.packages('stats')
	library(stats)
	if (!require(GenomicRanges)) install.packages('GenomicRanges')
	library(GenomicRanges)
	if (!require(igraph)) install.packages('igraph')
	library(igraph)
	if (!require(FNN)) install.packages('FNN')
	library(FNN)	
	if (!require(matrixcalc)) install.packages('matrixcalc')
	library(matrixcalc)	
	if (!require(matrixStats)) install.packages('matrixStats')
	library(matrixStats)	
	print("External packages loaded..")
} # end of external packages

###
print("Loading Paired-map functions...")
{ ### Paired-map classes
	{ ## distMatrix class
		methods::setClassUnion("MatrixOrmatrix", c("Matrix", "matrix"))
		distMatrix<-setClass(
			Class="distMatrix",
			slots=list(
				rmat.rna="MatrixOrmatrix",
				nmat.rna="MatrixOrmatrix",
				norm.rna="logical",
				method.rna="character",
				rmat.dna="MatrixOrmatrix",
				nmat.dna="MatrixOrmatrix",
				norm.dna="logical",
				method.dna="character",
				rmat.tag="MatrixOrmatrix",
				nmat.tag="MatrixOrmatrix",
				norm.tag="logical",
				method.tag="character",
				nmat.int="MatrixOrmatrix",
				int.method="character",
				int.input="character",
				p1="numeric",
				p2="numeric"
			)
		)
		setMethod(
			f='show',
			signature='distMatrix',
			definition=function(object){
				cat(
					"# of cells: ", nrow(object@rmat.rna),"\n",
					"# of RNA dims: ", ncol(object@rmat.rna), "\n",
					"Normalized RNA: ", object@norm.rna, "\n",
					"Method RNA: ", object@method.rna, "\n",
					"# of ATAC dims: ", ncol(object@rmat.dna), "\n",
					"Normalized DNA: ", object@norm.dna, "\n",
					"Method DNA: ", object@method.dna, "\n",
					"# of Tag dims: ", ifelse(is.null(ncol(object@rmat.tag)), "NA", ncol(object@rmat.tag)), "\n",
					"Normalized Tag: ", ifelse(is.null(object@norm.tag), "NA", object@norm.tag), "\n",
					"Method Tag: ", ifelse(is.null(object@method.tag), "NA", object@method.tag), "\n",
					"# of integrated dims: ", ifelse(is.null(object@imat), "NA", ncol(object@imat)), "\n",
					"Input of integrated: ", ifelse(is.null(object@imat.input), "NA", object@imat.input), "\n",
					"Method of integration: ", ifelse(is.null(object@imat.method), "NA", object@imat.method), "\n"
				)
			}
		)
		newDistMatrix<-function(){
			res=new("distMatrix",
				rmat.rna=matrix(0,0,0),
				nmat.rna=matrix(0,0,0),
				rmat.dna=matrix(0,0,0),
				nmat.dna=matrix(0,0,0),
				rmat.tag=matrix(0,0,0),
				nmat.tag=matrix(0,0,0),
				imat=matrix(0,0,0),
				imat.method=character(),
				imat.input=character(),
				method.rna=character(),
				method.dna=character(),
				method.tag=character(),
				p1=numeric(),
				p2=numeric(),
				norm.rna=FALSE,
				norm.dna=FALSE,
				norm.tag=FALSE
			)
		}
	} # end of distMatrix class
	{ ## dimReduction class
		dimReduction<-setClass(
			Class="dimReduction",
			slots=list(
				dmat.rna="matrix",
				sdev.rna="numeric",
				iter.rna="numeric",
				method.rna="character",
				dmat.dna="matrix",
				sdev.dna="numeric",
				iter.dna="numeric",
				method.dna="character",
				dmat.tag="matrix",
				sdev.tag="numeric",
				iter.tag="numeric",
				method.tag="character",
				dmat.int="matrix",
				sdev.int="numeric",
				iter.int="numeric",
				method.int="character"
				)
			)
		setMethod(
			f="show",signature="dimReduction",
			definition=function(object){
				cat(
					"RNA  - # of dims:", ncol(x=object@dmat.rna), "Method: ", object@method.rna, "\n",
					"ATAC - # of dims:", ncol(x=object@dmat.dna), "Method: ", object@method.dna, "\n",
					"Tag  - # of dims:", ncol(x=object@dmat.tag), "Method: ", object@method.tag, "\n"
					)}
			)
		newDimReduction<-function(object){
			res=new("dimReduction",
				dmat.rna=matrix(0,0,0),
				sdev.rna=numeric(),
				method.rna=character(),
				iter.rna=numeric(),
				dmat.dna=matrix(0,0,0),
				sdev.dna=numeric(),
				method.dna=character(),
				iter.dna=numeric(),
				dmat.tag=matrix(0,0,0),
				sdev.tag=numeric(),
				method.tag=character(),
				iter.tag=numeric(),
				dmat.int=matrix(0,0,0),
				sdev.int=numeric(),
				method.int=character(),
				iter.int=numeric()
				)
		}
	} # end of dimReduction class
	{ ## dataMatrix class
		dataMatrix<-setClass(
			Class="dataMatrix",
			slots=list(
				#mat.rna="dgTMatrix",
				mat.rna="Matrix",
				#genes="GRanges",
				genes="character",
				#mat.dna="dgTMatrix",
				mat.dna="Matrix",
				#bins.dna="GRanges",
				bins.dna="character",
				mat.tag="Matrix",
				#bins.tag="GRanges"
				bins.tag="character"
			)
		)
		setMethod(
			f="show",signature="dataMatrix",
			definition=function(object){
				cat(
					"# of cells:", ncol(x=object@mat.rna), "\n",
					"# of genes:", nrow(x=object@mat.rna), "\n",
					"# of ATAC bins:", nrow(x=object@mat.dna), "\n",
					"# of tag bins:", nrow(x=object@mat.tag), "\n",
					sep=""
				)
			}
		)
		newDataMatrix<-function(){
			mat.dna=matrix(0,0,0)
			#genes=GRanges()
			genes=as.character(c())
			mat.rna=matrix(0,0,0)
			#bins.dna=GRanges()
			bins.dna=as.character(c())
			mat.tag=matrix(0,0,0)
			#bins.tag=GRanges()
			bins.tag=as.character(c())
		}
	} # end of dataMatrix class
	{ ## embedding class
		methods::setClassUnion("MatrixOrmatrix", c("Matrix", "matrix"))
		embedding<-setClass(
			Class="embedding",
			slots=list(
				vis.dna="MatrixOrmatrix",
				vis.rna="MatrixOrmatrix",
				vis.tag="MatrixOrmatrix",
				vis.int="MatrixOrmatrix",
				method="character"
			)
		)
		setMethod(
			f="show",signature="embedding",
			definition=function(object){
				cat(
					"Method: ", object@method, "\n"
				)
			}
		)
		newEmbedding<-function(){
			vis.rna=matrix(nrow=0,ncol=0)
			vis.dna=matrix(nrow=0,ncol=0)
			vis.tag=matrix(nrow=0,ncol=0)
			vis.int=matrix(nrow=0,ncol=0)
			method=character()
		}
	} # end of embedding class
	{ ## cluster class 
		cluster<-setClass(
			Class="cluster",
			slots=list(
				id.rna="numeric",
				id.dna="numeric",	
				id.tag="numeric",
				id.int="numeric",
				method="character"
			)
		)
		setMethod(
			f="show",signature="cluster",
			definition=function(object){
				NULL ### need to be filled
			}
		)
		newCluster<-function(){
			id.rna=numeric()
			id.dna=numeric()
			id.tag=numeric()
			id.int=numeric()
			method=character()
		}
	} # end of cluster class
	{ ## paired class
		methods::setClassUnion("MatrixOrmatrix", c("Matrix", "matrix"))
		paired<-setClass(
			Class="paired",
			slots=list(
				project="character",	## project name
				barcode="character",	## cell id
				mat.filt="dataMatrix",## filtered data matrix	
				mat.raw="dataMatrix",	## joint data matrix	
				mat.nor="dataMatrix", ## normalized data matrix
				dis="distMatrix",			## Distance matrix 
				pca="dimReduction",		## dimension reduction
				tsne="embedding",			## cell embeddings
				umap="embedding",			## cell embeddings
				cluster="cluster",
				cell.downsample="numeric"
				)
			)
		setMethod("show", signature = "paired",
			definition=function(object){
				if((x=length(object@project))>=0L){
					cat("Paired-seq/tag dataset: ", object@project, "\n")
					#cat("# of cells: ", ifelse(is.null(length(object@barcode), 0, length(object@barcode))), "\n")
					cat("# of cells: ", length(object@barcode), "\n", sep="")
					cat("==Raw matrices==\n")
					cat("# of genes: ", nrow(object@mat.raw@mat.rna), "\n", sep="")
					cat("# of ATAC bins: ", nrow(object@mat.raw@mat.dna), "\n", sep="")
					cat("# of Tag bins: ", nrow(object@mat.raw@mat.tag), "\n", sep="")
					cat("==Filtered matrices==\n")
					cat("# of genes: ", nrow(object@mat.filt@mat.rna), "\n", sep="")
					cat("# of ATAC bins: ", nrow(object@mat.filt@mat.dna), "\n", sep="")
					cat("# of Tag bins: ", nrow(object@mat.filt@mat.tag), "\n", sep="")
					cat("==Normalized matrices==\n")
					cat("# of genes: ", nrow(object@mat.nor@mat.rna), "\n", sep="")
					cat("# of ATAC bins: ", nrow(object@mat.nor@mat.dna), "\n", sep="")
					cat("# of Tag bins: ", nrow(object@mat.nor@mat.tag), "\n", sep="")
				}
			}
		)
	} # end of paired class
} # end of Paired-map classes
###
{ ### Paired-map functions 
	{ ## newPaired 
		methods::setClassUnion("MatrixOrmatrix", c("Matrix", "matrix"))
		newPaired=function(){
			project=character()
			barcode=as.character(c())
			mat.filt=dataMatrix()
			mat.raw=dataMatrix()
			mat.nor=dataMatrix()
			dis=distMatrix()
			pca=dimReduction()
			tsne=embedding()
			umap=embedding()
			cluster=cluster()
			cell.downsample=as.numeric(c())
			res=new("paired",
			        project=project,
			        barcode=barcode,
			        mat.filt=mat.filt,
			        mat.raw=mat.raw,
			        mat.nor=mat.nor,
			        dis=dis,
			        pca=pca,
			        tsne=tsne,
			        umap=umap,
			        cluster=cluster,
			        cell.downsample=cell.downsample
			        )
		}
	} # end of newPaired

	{## loadPairedDataset
		loadPairedDataSet<-function(rna, dna, min.umi, min.genes, min.fragments, min.cells.rna, min.cells.dna, project){
			UseMethod("loadPairedDataSet")
		}
		loadPairedDataSet.default<-function(rna, dna, min.umi, min.genes, min.fragments, min.cells.rna, min.cells.dna, project){
			if(missing(rna) || missing(dna)){stop("Path to RNA and/or DNA is missing...")}
			if(missing(project)){project="Default Project"}
			if(missing(min.umi)){min.umi=200} # min transcript per nuclei 200
			if(missing(min.genes)){min.genes=100} # min genes per nuclei 100
			if(missing(min.fragments)){min.fragments=200} # min fragments per nuclei 200
			if(missing(min.cells.rna)){min.cells.rna=3} # min cells per gene 3
			if(missing(min.cells.dna)){min.cells.dna=3} # min cells per fragment 3
			print("Loading RNA and DNA matrices...")
			obj=newPaired()
			## Loading RNA matrix
			if(rna != "NA"){
				rna_mat=readMM(paste(rna, "matrix.mtx", sep="/"))
				rna_bc=read.csv(paste(rna,"barcodes.tsv",sep="/"), head=F)
				rna_genes=read.csv(paste(rna, "genes.tsv",sep="/"),sep="\t", head=F)
				colnames(rna_mat)<-rna_bc[,1]
				rownames(rna_mat)<-rna_genes[,2]
				rm(rna_bc)
				rm(rna_genes)
				## report RNA matrix stats
				print("Loaded RNA matrix")
				cat("# of barcodes: ", ncol(rna_mat), "\n", "# of genes: ", nrow(rna_mat), "\n\n", sep="")
				## pre filt RNA matrix
				n_umi<-colSums(rna_mat)
				rna_mat_filt<-rna_mat[,n_umi>min.umi]
				rm(rna_mat)
				rna_gene<-rna_mat_filt
				rna_gene@x[rna_gene@x>0]<-1
				n_gene<-colSums(rna_gene)
				n_cell_rna<-rowSums(rna_gene)
				rm(rna_gene)
				rna_mat_filt<-rna_mat_filt[n_cell_rna>min.cells.rna,n_gene>min.genes]
				cat("Min # of UMI: ", min.umi, "\n", "Min # of Genes: ", min.genes, "\n", "Min # of cells: ", min.cells.rna, "\n\n", sep="")
				cat("# of PF barcodes: ", ncol(rna_mat_filt), "\n", "# of PF genes: ", nrow(rna_mat_filt), "\n", "Median usable UMI: ", median(colSums(rna_mat_filt)), "\n\n", sep="")
				obj@mat.raw@genes=rownames(rna_mat_filt)
				obj@mat.raw@mat.rna=rna_mat_filt
				obj@barcode=colnames(rna_mat_filt)
			}


			## Loading DNA matrix
			if(dna != "NA"){
				dna_mat=readMM(paste(dna, "matrix.mtx", sep="/"))
				dna_bc=read.csv(paste(dna,"barcodes.tsv", sep="/"), head=F)
				dna_peaks=read.csv(paste(dna, "genes.tsv", sep="/"), head=F)
				colnames(dna_mat)<-dna_bc[,1]
				rownames(dna_mat)<-dna_peaks[,1]
				rm(dna_bc)
				rm(dna_peaks)
				## report DNA matrix ststs
				print("Loaded DNA matrix")
				cat("# of barcodes: ", ncol(dna_mat), "\n", "# of bins: ", nrow(dna_mat), "\n")
				## pre filt DNA matrix
				n_fragments<-colSums(dna_mat)
				dna_mat_filt<-dna_mat[,n_fragments>min.fragments]
				rm(dna_mat)
				dna_bins<-dna_mat_filt
				dna_bins@x[dna_bins@x>0]<-1
				n_cell_dna<-rowSums(dna_bins)
				rm(dna_bins)
				dna_mat_filt<-dna_mat_filt[n_cell_dna>min.cells.dna,]
				cat("Min # of fragments: ", min.fragments, "\n", "Min # of cells: ", min.cells.rna,"\n\n",sep="")
				cat("# of PF barcodes: ", ncol(dna_mat_filt), "\n", "# of PF bins: ", nrow(dna_mat_filt), "\n", "Median usable UMI: ", median(colSums(dna_mat_filt)),"\n\n",sep="")
				obj@mat.raw@bins.dna=rownames(dna_mat_filt)
				obj@mat.raw@mat.dna=dna_mat_filt
				obj@barcode=colnames(dna_mat_filt)
			}

			if(dna != "NA" & rna != "NA"){
				print("Reserving joint profiles...")
				joint_bc<-intersect(colnames(dna_mat_filt), colnames(rna_mat_filt))
				dna_mat_joint<-dna_mat_filt[,joint_bc]
				rna_mat_joint<-rna_mat_filt[,joint_bc]
				cat(
				    "# of joint cells: ", length(joint_bc), "\n",
				    "# of genes: ", nrow(rna_mat_joint), "\n",
				    "# of bins: ", nrow(dna_mat_joint), "\n",
				    "Median usable UMI: ", median(colSums(rna_mat_joint)), "\n",
				    "Median usable Fragments: ", median(colSums(dna_mat_joint)), "\n",
				    sep=""
				    )
				obj@mat.raw@genes=rownames(rna_mat_joint)
				obj@mat.raw@bins.dna=rownames(dna_mat_joint)
				colnames(rna_mat_joint)<-NULL
				rownames(rna_mat_joint)<-NULL
				colnames(dna_mat_joint)<-NULL
				rownames(dna_mat_joint)<-NULL
				obj@mat.raw@mat.rna=rna_mat_joint
				obj@mat.raw@mat.dna=dna_mat_joint
				obj@barcode=joint_bc
	
			}
			obj@project=project
			return(obj)
		}
		logMean<-function(input){
			mat.colsum<-base::colSums(input)
			mean.colsum<-mean(mat.colsum)
			scale.factor<-mat.colsum/mean.colsum
			mat.nor<-t(input)/scale.factor
			mat.nor<-log(mat.nor+1)
			mat.nor<-t(mat.nor)
			return(mat.nor)
		}
		logMedian<-function(input){
			mat.colsum<-colSums(input)
			mean.colsum<-median(mat.colsum)
			scale.factor<-mat.colsum/mean.colsum
			mat.nor<-t(input)/scale.factor
			mat.nor<-log(mat.nor+1)
			mat.nor<-t(mat.nor)
			return(mat.nor)
		}
		filtMatrix<-function(obj, input, high.threshold, low.threshold, method){
			UseMethod("filtMatrix", obj)
		}
		filtMatrix.default<-function(obj, input, high.threshold=5, low.threshold=-5, method){
			if(missing(obj)){stop("No Paired-seq dataset available...")}
			if(missing(input)){stop("Please define a matrix type.. dna/rna/tag")}
			if(missing(method)){stop("Please specify a method: binary/logMean")}
			input<-tolower(input)
			if(input=="dna"){raw_mat<-obj@mat.raw@mat.dna; features<-obj@mat.raw@bins.dna}
			if(input=="rna"){raw_mat<-obj@mat.raw@mat.rna; features<-obj@mat.raw@genes}
			if(input=="tag"){raw_mat<-obj@mat.raw@mat.rna; features<-obj@mat.raw@bins.tag}
			tmp_mat<-raw_mat
			tmp_mat@x[tmp_mat@x>0]<-1
			nCell_bin<-rowSums(tmp_mat)
			rm(tmp_mat)
			log_nCell_bin=log10(nCell_bin+1)
			s_log_nCell_bin=scale(log_nCell_bin)
			par(mar=c(2.5,3.5,3.5,2), lab=c(10,5,5))
			hist(s_log_nCell_bin, breaks=100, xlim=c(-5,5), col="grey", main="", xlab="Scale Features", tcl=-.25)
			val_bin<- c((s_log_nCell_bin>low.threshold)&(s_log_nCell_bin<high.threshold))
			filt_mat<-raw_mat[val_bin,]; filt_features<-features[val_bin]
			print("Filter features")
			cat("# total features: ", nrow(raw_mat), "\n", "# filtered features: ", nrow(filt_mat), "\n", sep="")
			cat("Normalizing...", " Method: ", method, "\n", sep="")
			if(method=="binary"){
				nor_mat<-filt_mat
				nor_mat@x[nor_mat@x>0]<-1
			}
			else if(method=="logMean"){
				nor_mat<-logMean(filt_mat)
			}else if(method=="none"){
				nor_mat<-filt_mat
			}else if(method=="logMedian"){
				nor_mat<-logMedian(filt_mat)
			}
			else{stop("Please specify a method: binary/logMean")}
			if(input=="dna"){obj@mat.filt@mat.dna<-filt_mat;obj@mat.filt@bins.dna=filt_features;obj@mat.nor@mat.dna<-nor_mat}
			if(input=="rna"){obj@mat.filt@mat.rna<-filt_mat;obj@mat.filt@genes=filt_features;obj@mat.nor@mat.rna<-nor_mat}
			if(input=="tag"){obj@mat.filt@mat.tag<-filt_mat;obj@mat.filt@bins.tag=filt_features;obj@mat.nor@mat.tag<-nor_mat}
			rm(raw_mat)
			rm(filt_mat)
			rm(filt_features)
			return(obj)
		}
	} # end of loadPairedSeq

	{	## calJaccard, adopted from SnapATAC
		calJaccard<-function(X_i, X_j){
			### input m = features
			### input n = cells
			### transpose here
			X_i<-t(X_i)
			X_j<-t(X_j)
			A = Matrix::tcrossprod(X_i, X_j)
			bi = Matrix::rowSums(X_i)
			bj = Matrix::rowSums(X_j)
			jmat = as.matrix(A / (replicate(ncol(A), bi) + t(replicate(nrow(A), bj)) - A))
			### output m = features
			### output n = cells
			### transpose here
			jmat<-t(jmat)
			return(jmat)				
		}
		calEuclidean<-function(mat.use){
			mat.use<-t(mat.use)
			smat<-apply(mat.use, 1, crossprod)
			mat1<-matrix(smat, nrow=dim(mat.use)[1], ncol=dim(mat.use)[1])
			mat3<-tcrossprod(mat.use)
			mat4<-mat1+t(mat1)-2*mat3
			diag(mat4)<-0
			mat5<-sqrt(mat4)
			mat5<-scale(mat5)
			return(mat5)
		}
		runDistance<-function(obj, input, feature.downsample, cell.downsample, seed.use, method){
			UseMethod("runDistance", obj)
		}
		runDistance.default<-function(obj, input, feature.downsample=1, cell.downsample=1, seed.use=131, method){
			if(missing(obj)){stop("No input object...")}
			if(missing(input)){stop("Please specify a input matrix...")};input<-tolower(input)
			if(missing(method)){stop("Please. specify a method... jaccard/euclidean")}
			if(length(obj@cell.downsample)>10L){
				cell.use<-obj@cell.downsample
				cat("Cell downsampling existed.\nUsing pre-defined cell subsampling...\n")
				not_new=1
			}
			else{
				if(cell.downsample>1 || cell.downsample<0){cell.downsample=1}
				if(cell.downsample==1){
					cell.use=c(1:length(obj@barcode))
				}
				else{
					cell.use.size=as.integer(length(obj@barcode)*cell.downsample)
					#cell.use<-c(rep(TRUE, cell.use.size), rep(FALSE, length(obj@barcode)-cell.use.size))
					set.seed(seed.use)
					cell.use=sample(length(obj@barcode), cell.use.size)
					#cell.use<-cell.use[sample(length(obj@barcode))]
					obj@cell.downsample=cell.use
				}
				not_new=0
			}
			if(feature.downsample>1 || feature.downsample<0){feature.downsample=1}
			if(input=="rna"){
				if(dim(obj@mat.nor@mat.rna)[1]<1L){
					stop("RNA matrix not normalized.")
				}
				mat.input=obj@mat.nor@mat.rna
			}
			if(input=="dna"){
				if(dim(obj@mat.nor@mat.dna)[1]<1L){
					stop("DNA matrix not normalized.")
				}
				mat.input=obj@mat.nor@mat.dna
			}
			if(input=="tag"){
				if(dim(obj@mat.nor@mat.tag)[1]<1L){
					stop("Tag matrix not normalized.")
				}
				mat.input=obj@mat.nor@mat.tag
			}
			mat.use<-mat.input
			if(feature.downsample == 1){
				mat.ref<-mat.use[,cell.use]
			}
			else{
				umi.mat.ref<-rowSums(mat.ref)
				feature.use.size=as.integer(feature.downsample*dim(mat.input)[1])
				if(length(umi.mat.ref[umi.mat.ref>0])<feature.ref.size){
					cat("Too few features with non-zero sum. Using all non-zero features.\n")
					feature.use.size=length(umi.mat.use[umi.mat.use>0])
				}
				feature.order.pool<-rev(order(umi.mat.ref))[1:feature.use.size]
				set.seed(seed.use)
				feature.use<-sample(feature.order.pool, feature.use.size)
				mat.use<-mat.input[feature.use,]
				mat.ref<-mat.input[feature.use,cell.use]
			}
			## check
			if(method=="jaccard"){
				cat("Calulating Jaccard Matrix...\n")
				jmat<-calJaccard(mat.use, mat.ref)
			}
			else if(method=="euclidean"){
				cat("Calulating Euclidean Matrix...\n")
				jmat<-calEuclidean(mat.use)
				jmat<-jmat[cell.use,]
				jmat<-as.matrix(jmat)
			}
			else{
				stop("Please. specify a method... jaccard/euclidean")
			}
			if(1){
				p1=Matrix::colMeans(mat.use)
				p2=p1[cell.use]
				obj@dis@p1=p1
				obj@dis@p2=p2				
			}
			if(input=="rna"){obj@dis@rmat.rna=jmat;obj@dis@method.rna=method}
			if(input=="dna"){obj@dis@rmat.dna=jmat;obj@dis@method.dna=method}
			if(input=="tag"){obj@dis@rmat.tag=jmat;obj@dis@method.tag=method}
			rm(jmat)
			rm(mat.use)
			rm(mat.ref)
			return(obj)
		}
	} # end of runDistance

	{	## runNormDistance
		.normOVE <- function(p1, p2){
		  pp = tcrossprod(p1, p2);
			ss = matrix(rep(p1,each=length(p2)), ncol=length(p2), byrow=TRUE) +  matrix(rep(p2, each=length(p1)), ncol=length(p2), byrow=FALSE)
			ee = pp/(ss - pp)
			ee<-t(ee)
			return(ee)	
		}
		trainRegressModel <- function(rmat, p1, p2){ ### adopted from SnapATAC
			idx.pairwise = which(rmat == 1, arr.ind=TRUE);
			emat = .normOVE(p1, p2);
			scale.factor = mean(rmat / emat);
			rmat[idx.pairwise] = scale.factor * emat[idx.pairwise];
			data = data.frame(x=c(emat), y=c(rmat));	
			model <- lm(y ~ x + I(x^2), data);
			return(model);	
		}
		normObservedJmat <- function(rmat, model, p1, p2, method){ ### adopted from SnapATAC

			idx.pairwise = which(rmat == 1, arr.ind=TRUE);
			emat = .normOVE(p1, p2);
			scale.factor = mean(rmat / emat);
			rmat[idx.pairwise] = scale.factor * emat[idx.pairwise];
			preds = predict(model, data.frame(x=c(emat)), se.fit = TRUE)
			if(method == "zscore"){
				norm = (c(rmat) - preds$fit) / (preds$se.fit);		
			}else if(method == "residual"){
				norm = c(rmat) -  preds$fit;
			}
			nmat = matrix(norm, nrow(emat), ncol(emat));	
			return(nmat);
		}

		runNormDistance<-function(obj, input,row.center, row.scale, low.threshold, high.threshold, seed.use, method){
			UseMethod("runNormDistance")
		}
		runNormDistance.default<-function(obj,input,row.center=TRUE,row.scale=TRUE,low.threshold=-5,high.threshold=5,seed.use=131,method="zscore"){
			if(missing(input)){stop("Please specify an input matrix")}
			if(input=="dna"){rmat=obj@dis@rmat.dna}
			if(input=="rna"){rmat=obj@dis@rmat.rna}
			if(input=="tag"){rmat=obj@dis@rmat.tag}
			p1=obj@dis@p1
			p2=obj@dis@p2

			if(method=="zscore" || method=="residual"){
				model.init=trainRegressModel(rmat, p1, p2) ### adopted from SnapATAC
				nmat=normObservedJmat(rmat, model.init, p1, p2, method=method) ### adopted from SnapATAC
			}
			if(method=="none"){
				nmat=rmat
			}
			if(row.center||row.scale){
				nmat=t(scale(t(nmat),center=row.center, scale=row.scale))
			}
			nmat[nmat>high.threshold]=high.threshold
			nmat[nmat<low.threshold]=low.threshold
			if(input=="dna"){obj@dis@nmat.dna=nmat}
			if(input=="rna"){obj@dis@nmat.rna=nmat}
			if(input=="tag"){obj@dis@nmat.tag=nmat}
			if(length(which(is.na(nmat)))>0){
				obj@cell.downsample<-c(0);
				cat("NaN lines appear, try to increase feature.downsample...\n")
				return(obj)
			}
			return(obj)
		}
	} # end of runNormDistance

	{ ## runFusionMatrices
		runFusionMatrices<-function(obj, use.rna, use.dna, use.tag, ratio,method){
			UseMethod("runFusionMatrices", obj)
		}
		runFusionMatrices.default<-function(obj, use.rna=T, use.dna=T, use.tag=F, ratio=c(1,1), method="hadamard"){
			if(missing(obj)){stop("Please specify an object...\n")}
			method=tolower(method)
			# specify input matrices
			use.three=FALSE
			if(use.rna && use.dna && !use.tag){
				mat1=minMaxScaling_non_zero(obj@dis@nmat.rna)
				if(obj@dis@method.rna=="jaccard"){mat1=1-mat1}
				mat2=minMaxScaling_non_zero(obj@dis@nmat.dna)
				if(obj@dis@method.dna=="jaccard"){mat2=1-mat2}
				cat("Fuse RNA and DNA matrices...\n")
				obj@dis@int.input="RNA+DNA"
			}else if(use.rna && !use.dna && use.tag){
				mat1=minMaxScaling_non_zero(obj@dis@nmat.rna)
				if(obj@dis@method.rna=="jaccard"){mat1=1-mat1}
				mat2=minMaxScaling_non_zero(obj@dis@nmat.tag)
				if(obj@dis@method.tag=="jaccard"){mat2=1-mat2}
				cat("Fuse RNA and Tag matrices...\n")
				obj@dis@int.input="RNA+Tag"
			}else if(!use.rna && use.dna && use.tag){
				mat1=minMaxScaling_non_zero(obj@dis@nmat.dna)
				if(obj@dis@method.dna=="jaccard"){mat1=1-mat1}
				mat2=minMaxScaling_non_zero(obj@dis@nmat.tag)
				if(obj@dis@method.tag=="jaccard"){mat2=1-mat2}
				cat("Fuse DNA and Tag matrices...\n")
				obj@dis@int.input="DNA+Tag"
			}else if(use.rna && use.dna && use.tag){
				mat1=minMaxScaling_non_zero(obj@dis@nmat.rna)
				if(obj@dis@method.rna=="jaccard"){mat1=1-mat1}
				mat2=minMaxScaling_non_zero(obj@dis@nmat.dna)
				if(obj@dis@method.dna=="jaccard"){mat2=1-mat2}
				mat3=minMaxScaling_non_zero(obj@dis@nmat.tag)
				if(Obj@dis@method.tag=="jaccard"){mat3=1-mat3}
				use.three<-TRUE
				cat("Fuse three matrices...\n")
				obj@dis@int.input="RNA+DNA+Tag"
			}
			if(use.three){mat3<-minMaxScaling_non_zero(mat3)}
			if(method=="hadamard"){
				nmat.int<-hadamard.prod(mat1, mat2)
				if(use.three){nmat.int<-hadamard.prod(nmat.int, mat3)}
				obj@dis@int.method="hadamard"
			}
			else if(method=="sum"){
				nmat.int<-mat1*ratio[1]+mat2*ratio[2]
				if(use.three){nmat.int<-nmat.int+mat3*ratio[3]}
				obj@dis@int.method="sum"
			}
			else{
				stop(method, "not supported..", sep=" ")
			}
			obj@dis@nmat.int=nmat.int
			return(obj)
		}
	} # end of runFusionMatrices

	{ ## pca
		runPCA<-function(obj,input,n){
			UseMethod("runPCA")
		}
		runPCA.default<-function(obj,input,n=50){
			if(missing(obj)){stop("Are you kidding me?")}
			if(missing(input)){stop("Please specify an input matrix")};input<-tolower(input)
			if(input=="dna"){mat=obj@dis@nmat.dna}
			if(input=="rna"){mat=obj@dis@nmat.rna}
			if(input=="tag"){mat=obj@dis@nmat.tag}
			if(input=="int"){mat=obj@dis@nmat.int}
			pca<-prcomp_irlba(t(mat), n=n)
			if(input=="dna"){
				obj@pca@dmat.dna<-pca$x
				obj@pca@sdev.dna<-pca$sdev
				obj@pca@iter.dna<-0
				obj@pca@method.dna<-"PCA"
			}
			if(input=="rna"){
				obj@pca@dmat.rna<-pca$x
				obj@pca@sdev.rna<-pca$sdev
				obj@pca@iter.rna<-0
				obj@pca@method.rna<-"PCA"
			}
			if(input=="tag"){
				obj@pca@dmat.tag<-pca$x
				obj@pca@sdev.tag<-pca$sdev
				obj@pca@iter.tag<-0
				obj@pca@method.tag<-"PCA"
			}
			if(input=="int"){
				obj@pca@dmat.int<-pca$x
				obj@pca@sdev.int<-pca$sdev
				obj@pca@iter.int<-0
				obj@pca@method.int<-"PCA"
			}
			return(obj)
		}
	}

	{ ## runUMAP
		runUMAP<-function(obj, input, k, use.dims,...){
			UseMethod("runUMAP")
		}
		runUMAP.default<-function(obj, input, k=15, use.dims=c(1,10),...){
			if(missing(obj)){stop("Are you kidding me?")}
			if(missing(input)){stop("Please specify an input matrix")};input<-tolower(input)
			if(input=="dna"){mat=obj@pca@dmat.dna[,use.dims]}
			if(input=="rna"){mat=obj@pca@dmat.rna[,use.dims]}
			if(input=="tag"){mat=obj@pca@dmat.tag[,use.dims]}
			if(input=="int"){mat=obj@pca@dmat.int[,use.dims]}
			umap.out<-umap(mat, n_neighbors=k, verbose=T, n_threads=8, ...)
			if(input=="dna"){obj@umap@vis.dna<-umap.out}
			if(input=="rna"){obj@umap@vis.rna<-umap.out}
			if(input=="tag"){obj@umap@vis.tag<-umap.out}
			if(input=="int"){obj@umap@vis.int<-umap.out}
			return(obj)
		}
		runTSNE<-function(obj, input, iter, use.dims,...){
			UseMethod("runTSNE")
		}
		runTSNE.default<-function(obj, input, iter=500, use.dims=c(1,10),...){
			if(missing(obj)){stop("Are you kidding me?")}
			if(missing(input)){stop("Please specify an input matrix")};input<-tolower(input)
			if(input=="dna"){mat=obj@pca@dmat.dna[,use.dims]}
			if(input=="rna"){mat=obj@pca@dmat.rna[,use.dims]}
			if(input=="tag"){mat=obj@pca@dmat.tag[,use.dims]}
			if(input=="int"){mat=obj@pca@dmat.int[,use.dims]}
			tsne.out<-Rtsne(mat, n_neighbors=k, pca=FALSE, verbose=T, num_threads=8, max_iter=iter, ...)
			if(input=="dna"){obj@tsne@vis.dna<-tsne.out$Y}
			if(input=="rna"){obj@tsne@vis.rna<-tsne.out$Y}
			if(input=="tag"){obj@tsne@vis.tag<-tsne.out$Y}
			if(input=="int"){obj@tsne@vis.int<-tsne.out$Y}
			return(obj)
		}

	}

	{ ## cluster
		runCluster<-function(obj,input,seed.use,use.dims,k,...){
			UseMethod("runCluster", obj)
		}
		runCluster.default<-function(obj,input,seed.use=131,use.dims=c(1,10),k=15,...){
			if(missing(obj)){stop("No object...\n")}
			if(missing(input)){stop("Please specify an input matrix")};input<-tolower(input)
			if(input=="dna"){mat=obj@pca@dmat.dna[,use.dims]}
			if(input=="rna"){mat=obj@pca@dmat.rna[,use.dims]}
			if(input=="tag"){mat=obj@pca@dmat.tag[,use.dims]}
			if(input=="int"){mat=obj@pca@dmat.int[,use.dims]}
			knn.norm=get.knn(as.matrix(mat), k=k)
			knn.norm = data.frame(from = rep(1:nrow(knn.norm$nn.index),  k), to = as.vector(knn.norm$nn.index), weight = 1/(1 + as.vector(knn.norm$nn.dist)))
			nw.norm = graph_from_data_frame(knn.norm, directed = FALSE)
			nw.norm = simplify(nw.norm)
			lc.norm.combine = cluster_louvain(nw.norm)
			if(input=="dna"){obj@cluster@id.dna=(lc.norm.combine$membership)}
			if(input=="rna"){obj@cluster@id.rna=(lc.norm.combine$membership)}
			if(input=="tag"){obj@cluster@id.tag=(lc.norm.combine$membership)}
			if(input=="int"){obj@cluster@id.int=(lc.norm.combine$membership)}
			return(obj)
		}
	} # end of cluster
	col.geneExpr<-colorRampPalette(c("gray85", "forestgreen", "darkgreen"))(10)
	col.promAccs<-colorRampPalette(c("gray85", "deepskyblue3", "deepskyblue4"))(10)
	col.depth<-colorRampPalette(c("grey85", "salmon1", "salmon4"))(10)
	t_col <- function(color, percent = 50, name = NULL) {
		rgb.val <- col2rgb(color)
		t.col <- rgb(rgb.val[1,], rgb.val[2,], rgb.val[3,],
             max = 255,
             alpha = (100 - percent) * 255 / 100,
             names = name)
		return(t.col)
	}
	minMaxScaling_non_zero<-function(x){
		return((x+0.001-min(x))/(max(x)+0.001-min(x)))
	}

	minMaxScaling<-function(x){
		return((x-min(x))/(max(x)-min(x)))
	}

	{ ## plotFeature
		plotFeature<-function(obj, feature, feature.type, norm.frac, norm.log,embedding.use, embedding.type,outlies, title,pch,cex, legend,legend.position,...){
			UseMethod("plotFeature", obj)
		}
		plotFeature.default<-function(obj, feature, feature.type="rna", norm.frac=FALSE, norm.log=FALSE, outliers=c(0.01,0.99),embedding.use="umap",pch=19, cex=0.25, title,embedding.type=feature.type,legend.position="right",legend=F,...){
			if(missing(obj)){stop("Please specify an object..\n")}
			if(missing(feature)){stop("Please specify an feature...\n")}
			feature.type<-tolower(feature.type)
			if(feature=="cluster"){
				if(feature.type=="dna"){
					feature.value=obj@cluster@id.dna
				}else if(feature.type=="rna"){
					feature.value=obj@cluster@id.rna
				}else if(feature.type=="tag"){
					feature.value=obj@cluster@id.tag
				}else if(feature.type=="int"){
					feature.value=obj@cluster@id.int
				}
				col=colPanel
				legend=T
			}else	if(feature=="depth"){
				if(feature.type=="dna"){
					feature.value=base::colSums(obj@mat.raw@mat.dna)
				}else if(feature.type=="rna"){
					feature.value=base::colSums(obj@mat.raw@mat.rna)
				}else if(feature.type=="tag"){
					feature.value=base::colSums(obj@mat.raw@mat.tag)
				}
				feature.value=log(feature.value+1)
				col<-col.depth
				feature.value<-minMaxScaling(feature.value)*9+1
			}else{
				### plot gene expression or promoter accessibility
				if((feature.type)=="rna"){
					idx=which(obj@mat.raw@genes==feature)
					if(length(idx)==0){stop("Cannot find ", feature, "\n", sep="")}
					feature.value=obj@mat.raw@mat.rna[idx,]
					col=col.geneExpr
				}else	if((feature.type)=="dna"){
					idx=which(obj@mat.raw@bins.dna==feature)
					if(length(idx)==0){stop("Cannot find ", feature, "\n", sep="")}
					feature.value=obj@mat.raw@mat.rna[idx,]
					col=col.promAccs
				}
				if(norm.frac){
					if(feature.type=="rna"){
						cell.size<-base::colSums(obj@mat.raw@mat.rna)
					}else if(feature.type=="dna"){
						cell.size<-base::colSums(obj@mat.raw@mat.dna)
					}
					feature.value<-feature.value/cell.size
				}
				## log normalize
				if(norm.log){
					feature.value=log(feature.value+1)
				}
								## filter outliers
				outlier.low=quantile(feature.value, outliers[1])
				outlier.high=quantile(feature.value, outliers[2])
				feature.value[feature.value>outlier.high]=outlier.high
				feature.value[feature.value<outlier.low]=outlier.low
				feature.value<-minMaxScaling(feature.value)*9+1
			}
			if(missing(title)){title=paste(toupper(feature.type), feature, "on", toupper(embedding.type), toupper(embedding.use), sep=" ")}
			if(tolower(embedding.use)=="umap"&&tolower(embedding.type)=="rna"){mat=obj@umap@vis.rna;labx="UMAP1";laby="UMAP2"}
			if(tolower(embedding.use)=="umap"&&tolower(embedding.type)=="dna"){mat=obj@umap@vis.dna;labx="UMAP1";laby="UMAP2"}
			if(tolower(embedding.use)=="umap"&&tolower(embedding.type)=="tag"){mat=obj@umap@vis.tag;labx="UMAP1";laby="UMAP2"}
			if(tolower(embedding.use)=="umap"&&tolower(embedding.type)=="int"){mat=obj@umap@vis.int;labx="UMAP1";laby="UMAP2"}
			if(tolower(embedding.use)=="tsne"&&tolower(embedding.type)=="rna"){mat=obj@tsne@vis.rna;labx="tSNE1";laby="tSNE2"}
			if(tolower(embedding.use)=="tsne"&&tolower(embedding.type)=="dna"){mat=obj@tsne@vis.dna;labx="tSNE1";laby="tSNE2"}
			if(tolower(embedding.use)=="tsne"&&tolower(embedding.type)=="tag"){mat=obj@tsne@vis.tag;labx="tSNE1";laby="tSNE2"}
			if(tolower(embedding.use)=="tsne"&&tolower(embedding.type)=="int"){mat=obj@tsne@vis.int;labx="tSNE1";laby="tSNE2"}
			ncol=0
			if(legend){ncol=as.integer(max(feature.value/20)+1)} 
			par(mar=c(2.5, 2.5, 3, 3.0+1.20*ncol),xpd=TRUE, col.axis="black", cex.axis=1, lab=c(3,3,5), mgp=c(1.25,0.25,0))
			plot(mat, pch=pch,cex=cex,col=col[feature.value], main=title,xlab=labx, ylab=laby, bty="l", tcl=-.25, xgap.axis=2, ygap.axis=2, font.lab=2, ...)
			if(legend){legend("right", inset=c(-0.12,0),legend=c(1:max(feature.value)), pch=19, cex=0.75, ncol=ncol, bty="n", col=col)}
		}

		plotPCA<-function(obj, input){
			UseMethod("plotPCA", obj)
		}
		plotPCA.default<-function(obj, input="dna"){
			#opar<-par(newsettings)
			par(mar=c(2.5,3.5,3.5,2), lab=c(10,5,5))
			if(missing(obj)){stop("Please specify an object..")}
			layout(matrix(c(1,1,2,3,4,1,1,5,6,7,8,9,10,11,12,13,14,15,16,17),4,5,byrow=TRUE))
			if(input=="dna"){pca<-obj@pca@dmat.dna;sdev<-obj@pca@sdev.dna}	
			if(input=="rna"){pca<-obj@pca@dmat.rna;sdev<-obj@pca@sdev.rna}	
			if(input=="int"){pca<-obj@pca@dmat.int;sdev<-obj@pca@sdev.int}	
			if(input=="tag"){pca<-obj@pca@dmat.tag;sdev<-obj@pca@sdev.tag}	
			
			plot(sdev, pch=19, col="grey", cex=0.5, xlab="PC", ylab="sdev", main=paste("PCA sdev", sep=" "), tcl=-.25)
			barcodes<-obj@barcode
			write.table(barcodes, file="bc.tmp", quote=F, col.names=F, row.names=F)
			bc.split<-read.csv("bc.tmp", sep=":", head=F)
			feature.value<-bc.split[,dim(bc.split)[2]]
			col<-colPanel
			for(i in 1:16){
				pc1<-i*2-1;pc2<-i*2
				plot(pca[,pc1:pc2], pch=19, cex=0.05, col=colPanel[feature.value], tcl=-.25, main=paste("PC", pc1, " vs PC", pc2, sep=""))
			}
			par(mfrow=c(1,1))
		}
	} # end of plotFeature

	{ ## plotPairedSeq 
	} # end of plotPairedSeq
	{ ## plotPairedTag 
	} # end of plotPairedTag	

	colPanel = c(
		"#E31A1C", "#FFD700", "#771122", "#777711", "#1F78B4", "#68228B", "#AAAA44",
		"#60CC52", "#771155", "#DDDD77", "#774411", "#AA7744", "#AA4455", "#117744", 
		"#000080", "#44AA77", "#AA4488", "#DDAA77", "#D9D9D9", "#BC80BD", "#FFED6F",
	    "#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#BF5B17",
	    "#666666", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02",
	    "#A6761D", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C",
	    "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#B15928", "#FBB4AE", "#B3CDE3",
	    "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC", "#F2F2F2",
	    "#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4", "#E6F5C9", "#FFF2AE", "#F1E2CC",
	    "#CCCCCC", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FFFF33", "#A65628",
	    "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854",
	    "#FFD92F", "#E5C494", "#B3B3B3", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072",
	    "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5"
  )


} # end of Paired-map functions





print("Finished loading Paired-map functions...")
print(paste("Version: ", version, sep=""))















































































