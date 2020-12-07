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

*** Please modification the paths to reference files according to the annotations in the script file.

#### The output of this step includes:

<code>Sample_combined.fq.gz</code>  This file is a combined fastq file including Read1 sequences/qualities and barcode sequences extracted from Read2.

<code>Sample_BC.sam</code> This is a temporally file used to assign extracted barcode sequences to Cellular Barcode. Please delete this file if you have successful obtain <code>Sample_BC_cov.fq.gz</code>.

<code>Sample_BC_cov.fq.gz</code> This is the fastq file with Read1 sequences and qualities, the Cellular Barocde and UMI from Read2 are now in ReadName section of the fastq file (and subsequent alignment files).

The Cellular Barcode and UMI are in the format of <code>ILLUMINA_READ_NAME:aa:bb:cc:UMI</code>

<code>aa</code> and <code>bb</code> are numeric values between <code>01</code> and <code>96</code> representing the wells of the 2nd and 1st round of barcode ligation. <code>cc</code> is a numeric value between <code>01</code> and <code>12</code> indicating the 12 tubes for multiplexing samples during tagmentation a nd reverse transcription steps.

#### The following metrics from this stetp can be used for QC:

Output from <code>reachtools combine2</code> step:

<code>205158	229195	89.51% of reads have full barcodes for Test sample.</code>. 

This ratio if the percentage of reads that can sucessfully extract all 3 barcodes from Read2.

Typically, >85% and >75% of reads from DNA and RNA libraries will have all 3 barcodes.

Output from <code>reachtools convert2</code> step:

<code>205158 reads processed.</code>

<code>185824 mapped reads.</code>

This numbers are the reads that can be uniquely assigned to one Cellular Barocdes.

Typically, >85% of reads (both DNA and RNA) can be assigned.


## 2 Mapping to the genome
For DNA reads, we used bowtie2; for RNA reads, we used STAR.

Use shellscrips/02.proc_DNA.sh and shellscrips/03.proc_RNA.sh.

After these, the fastq files of each sub-libraries were then converted to cell-counts matrices.

*** Please modifiy the variables in the scripts according to your environment.

Typically, >85% of RNA reads can be mapped to reference genome using STAR, but >60% is also accepatble. For different histone marks, from 60% - 95% of DNA reads can be mapped to reference genome using bowtie2 (including reads mapped to multiple loci).


## 3 Merge sub-libraries for downstream analysis
The last round of combinatorial index is PCR indexing (sub-libraries). 

To merge matrices from different sub-libraires, an unique prefix should added to the cellular barcodes for each sub-library. Please also make sure the DNA and RNA sublibraries share the same set of sub-library-specific prefixes.

Use perlscripts/merge_mtx.pl. Please see annotations in the script file for details.

## 4 Downstream custom analyses
You can cluster the single cells for DNA and RNA independely or jointly. 

For RNA analysis, Seurat (http://www.satijalab.org/seurat) is recommended. 

For DNA analysis, snapATAC (https://github.com/r3fang/SnapATAC) is recommeded.

You can also use Paired-map (https://github.com/cxzhu/Paired-map) to do the clustering. Paired-map contain few essential functions for analysis of Paired-seq/tag dataset but still under development. The Jaccard module was borrowed from snapATAC.
