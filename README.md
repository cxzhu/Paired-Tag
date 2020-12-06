# Paired-Tag Analysis
#### Please have the following softwares installed first:

bowtie, http://bowtie-bio.sourceforge.net/index.shtml

bowtie2, http://bowtie-bio.sourceforge.net/bowtie2/index.shtml

Trim_galore, https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/

STAR, https://github.com/alexdobin/STAR


#### Please follow the following steps to prepare the scripts for Paired-Tag data preprocessing.

1. Download and uncompress the scripts: 

<code>wget https://github.com/cxzhu/Paired-Tag/archive/master.zip</code>

<code>unzip master.zip</code>

2. Build the "reachtools" tool:

<code>cd reachtools</code>

<code>sh make.sh</code>

3. Build the Cellular Barocdes reference:

<code>cd ../references</code>

<code>bowtie-build ./cell_id_full.fa cell_id</code>



#### Analysis of Paired-Tag/Paired-seq datasets include the following steps:


## 1 Pre-processing
Extract cellular barcode from Read2, map the reads to reference cell_ID, and convert the mapped cell ID samfiles to useable fastq files.

Use shellscrips/01.pre_process_paired_tag_fastq.sh.

*** See the comments in this script if you are processing fastq files downloaded directly from GEO.

## 2 Mapping to the genome
For DNA reads, we used bowtie2; for RNA reads, we used STAR.

Use shellscrips/02.proc_DNA.sh and shellscrips/03.proc_RNA.sh.

After these, the fastq files of each sub-libraries were then converted to cell-counts matrices.


## 3 Merge sub-libraries for downstream analysis
The last round of combinatorial index is PCR indexing (sub-libraries). 

To merge matrices from different sub-libraires, an unique prefix should added to the cellular barcodes for each sub-library. Please also make sure the DNA and RNA sublibraries share the same set of sub-library-specific prefixes.

Use perlscripts/merge_mtx.pl.

## 4 Downstream custom analyses
You can cluster the single cells for DNA and RNA independely or jointly. 

For RNA analysis, Seurat (http://www.satijalab.org/seurat) is recommended. 

For DNA analysis, snapATAC (https://github.com/r3fang/SnapATAC) is recommeded.

You can also use Paired-map (https://github.com/cxzhu/Paired-map) to do the clustering. Paired-map contain few essential functions for analysis of Paired-seq/tag dataset but still under development. The Jaccard module was borrowed from snapATAC.
