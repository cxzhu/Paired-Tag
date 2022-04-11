# Preprocessing of Paired-Tag/Paired-seq datasets

#### For latest protocol, please refer to [Paired-Tag protocol & FAQs](https://github.com/cxzhu/Paired-Tag/tree/master/protocol).

#### Please have the following softwares installed first:

- bowtie, http://bowtie-bio.sourceforge.net/index.shtml
   
   For bowtie version=1.x, please modify <code>shellscrips/01.pre_process_paired_tag_fastq.sh</code> according to the comments in that file.

- bowtie2, http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
- samtools, http://www.htslib.org/
   samtools version >= 1.3.1 is required.

- Trim_galore, https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/

- STAR, https://github.com/alexdobin/STAR

- Optional: FastQC, https://www.bioinformatics.babraham.ac.uk/projects/fastqc/


#### Initial QC
A quick initial QC can be done will FastQC software. FastQC will give a QC summary for several key quality metrics of fastq files generated from Illumina bcl2fastq program. The image below shows the "Per base sequence content" and "Adapter Concent" sections of the FastQC output file from a representative Paired-Tag library. 

As shown in Read1 report, there is a high fraction of G base in the 2nd base, that is expected from the library construction (no trimming are needed for this part as Bowtie2 will handle it properly). For read2, there are 3 base balanced regions (UMI and barcode) between 3 base-inbalenced linkers (as indicated in the image). If the linker regions does not show high fluctuation as in the representative image, that may indicate a low ligation efficiency.

For "Adapter Content" section, typically we will expect a low fraction of Nextera adaptor sequence (at 100th bp, 5%-20%; expect higher percentage if sequenced to 150 bp or longer) in Read2 library and negligible adaptor content from Read1 library. Higher fraction of adaptors (1) in RNA library indicates: amount of N5-Tn5 is too high in 2nd adaptor tagging step of RNA library, or (2) in DNA library: tagmentation efficiency (antibody efficiency) are low, may expect low library complexity.

![Image_of_QC](https://github.com/cxzhu/Paired-Tag/blob/master/img/QC.png)

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
or
<code>bowtie-build ./cell_id_full_407.fa cell_id_407</code>

The two references provided are different in the barcodes length. In order to increase the throughput of Paired-Tag, we further scaled up the number of barcodes used for 2nd and 3rd round barocding from 96 to 384. To ensure hamming distance between barcodes we therefore increase barcodes length from 7bp to 8bp. If you are using your own barocdes please genreate the barcodes reference based on your sequence. 


#### Analysis of Paired-Tag/Paired-seq datasets include the following steps:


## 1. Pre-processing
Extract cellular barcode from Read2, map the reads to reference cell_ID, and convert the mapped cell ID samfiles to useable fastq files.

*** Please modification the paths to reference files according to the annotations in the script file.

Use <code>shellscrips/01.pre_process_paired_tag_fastq.sh</code>.

*** As SRA will trim the header of fastq files and thus a different script will be needed. See the annotations in this file if you are processing fastq files downloaded from [SRA](https://www.ncbi.nlm.nih.gov/sra?term=SRP266461).


#### The output of this step includes:

1. <code>Sample_combined.fq.gz</code>  This file is a combined fastq file including Read1 sequences/qualities and barcode sequences extracted from Read2.

2. <code>Sample_BC.sam</code> This is a temporally file used to assign extracted barcode sequences to Cellular Barcode. Please delete this file if you have successful obtained <code>Sample_BC_cov.fq.gz</code>.

3. <code>Sample_BC_cov.fq.gz</code> This is the fastq file with Read1 sequences and qualities, the Cellular Barocde and UMI from Read2 are now in ReadName section of the fastq file (and subsequent alignment files).

The Cellular Barcode and UMI are in the format of <code>ILLUMINA_READ_NAME:aa:bb:cc:UMI</code>

<code>aa</code> and <code>bb</code> are representing the wells of the 2nd and 1st round of barcode ligation. <code>cc</code> is a numeric value between <code>01</code> and <code>12</code> indicating the 12 tubes for multiplexing samples during tagmentation a nd reverse transcription steps.

#### The following metrics from this stetp can be used for QC:

Report from <code>reachtools combine2</code> step:

<code>205158	229195	89.51% of reads have full barcodes for Test sample.</code>. 

This ratio if the percentage of reads that can sucessfully extract all 3 barcodes from Read2.

Typically, >85% and >75% of reads from DNA and RNA libraries will have all 3 barcodes.

Report from <code>reachtools convert2</code> step:

<code>205158 reads processed.</code>

<code>185824 mapped reads.</code>

This numbers are the reads that can be uniquely assigned to one Cellular Barocdes.

Typically, >85% of reads (both DNA and RNA) can be assigned.


## 2. Mapping to the genome
For DNA reads, we used bowtie2; for RNA reads, we used STAR.

Use <code>shellscrips/02.proc_DNA.sh</code> and <code>shellscrips/03.proc_RNA.sh</code>.

After these, the fastq files of each sub-libraries were then converted to cell-counts matrices.

*** Please modifiy the variables in the scripts according to your environment.

Typically, >85% of RNA reads can be mapped to reference genome using STAR, but >60% is also accepatble. For different histone marks, from 60% - 95% of DNA reads can be mapped to reference genome using bowtie2 (including reads mapped to multiple loci).


## 3. Merge sub-libraries for downstream analysis
The last round of combinatorial index is PCR indexing (sub-libraries). 

Optional but recommended: Filtering barcode with low reads numbers.
Before merging sub-libraries, it is recommended to count the DNA and RNA reads numbers for each sub-libraries seperatedly and filter the low quality barcodes. A simple way is to plot the # of DNA and RNA reads for each barcodes (x-y in log scale, as in the picture shown below) and find a suitable cutoff from the plot, use <code>rscripts/plot_reads_numbers.R</code>.

![image_of_reads](https://github.com/cxzhu/Paired-Tag/blob/master/img/reads_plot.png)

Save the ID of barcodes with both high DNA and RNA reads number in a seperate file, and the <code>perlscript/filt_mtx.pl</code> to filter the raw matrix.

To merge matrices from different sub-libraires, an unique prefix should added to the cellular barcodes for each sub-library. Please also make sure the DNA and RNA sublibraries share the same set of sub-library-specific prefixes.

Use <code>perlscripts/merge_mtx.pl</code>. Please see annotations in the script file for details.

After merging, the sub-library ID will be added into the front of Cellular Barcodes, turns to <code>ii:aa:bb:cc</code>, where <code>ii</code> is the sub-library ID defined in <code>merge_list.txt</code> (see the picture below). Please use the same set of sub-library IDs for the paired DNA and RNA datasets.

![Image_of_Barcode_format](https://github.com/cxzhu/Paired-Tag/blob/master/img/barcodes_format-01.png)

## 4. Downstream custom analyses
You can cluster the single cells for DNA and RNA independely or jointly. The processed matrices files are available from Supplementary file of this [GEO series](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152020). E.g., for RNA analysis, Seurat (http://www.satijalab.org/seurat) can be used for scRNA-seq and  snapATAC (https://github.com/r3fang/SnapATAC) can used for DNA analysis.


<br>
<br>
  
For more information, please refer to:

- [Joint profiling of histone modifications and transcriptome in single cells from mouse brain.](https://www.nature.com/articles/s41592-021-01060-3) <em>Nature Methods</em>, 2021

A step-by-step protocol of Paired-Tag is also available from:

- [High-throughput single-cell joint analysis of histone modifications and gene expression by Paired-Tag.](https://protocolexchange.researchsquare.com/article/pex-1301/v1)



