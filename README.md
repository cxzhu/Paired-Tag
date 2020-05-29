# Paired-tag Analysis
Analysis of Paired-tag datasets include the following steps:

## 1 Pre-processing
Extract cellular barcode from Read2, map the reads to reference cell_ID, and convert the mapped cell ID samfiles to useable fastq files.

Use shellscrips/01.pre_process_paired_tag_fastq.sh.

## 2 Mapping to the genome
For DNA reads, we used bowtie2; for RNA reads, we used STAR.

Use 
shellscrips/02.proc_DNA.sh
and 
shellscrips/03.proc_RNA.sh
.

After these, the fastq files of each sub-libraries were then converted to cell-counts matrices.

## 3 Merge sub-librareis for down-stream analysis
The last round of combinatorial index is PCR indexing (sub-libraries). 

To merge matrices from different sub-libraires, an unique prefix should added to the cellular barcodes for each sub-library. Please also make sure the DNA and RNA sublibraries share the same set of sub-library-specific prefixes.

Use perlscripts/merge_mtx.pl.

## 4 Downstream custom analyses
You can cluster the single cells for DNA and RNA independely or jointly. 

For RNA analysis, Seurat (http://www.satijalab.org/seurat) is recommended. 

For DNA analysis, snapATAC (https://github.com/r3fang/SnapATAC) is recommeded.

You can also use Paired-map (https://github.com/cxzhu/Paired-map) to do the clustering. Paired-map contain few essential functions for analysis of Paired-seq/tag dataset but still under development. The Jaccard module was borrowed from snapATAC.
