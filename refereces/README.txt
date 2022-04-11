This folder includes the references used for analysis of mm10 and hg38 data. For other version of reference genome, please make your own ref file with the same format.

- cell_id_full.fa ## used to generate a 7bp-length, 96 barcodes bowtie index for demultiplexing cellular barcodes
- cell_id_full_407.fa ## used to generate a 8bp-length, 384 barcodes bowtie index for demultiplexing cellular barcodes. Barcodes length is also different from cell_id_full.fa 

- mm10.RNA.txt ## mouse ref file for generating cell-to-genes matrix with `reachtools bam2Mtx2` function
- mm10.bin5k.txt ## mouse ref file for generating cell-to-bins matrix with `reachtools bam2Mtx2` function
- hg38.RNA.txt ## human ref file (GENCODE v38) for generating cell-to-genes matrix with `reachtools bam2Mtx2` function
- hg38.bin5k.txt ## human ref file for generating cell-to-bins matrix with `reachtools bam2Mtx2` function
