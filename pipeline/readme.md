## Steps

### 1 Customize environment
a Install all packages needed, including bowtie, bowtie2, STAR, trim_galore, samtools

b Download all references and scripts from https://github.com/cxzhu/Paired-Tag

c Modify the runtime path in "run.sh"

### 2 Run pipeline
a Modify the sample information in filehead of "run.sh"

b Rename or linke the fastq files to "${prefix}${ID}_R1.fq.gz" and "${prefix}${ID}_R2.fq.gz", or "${prefix}${ID}_R1.fastq.bz2" and "${prefix}${ID}_R2.fastq.bz2"

c Run preprocessing by typing "sh run.sh Pre_bz2" (for *fastq.bz2 files) or "sh run.sh Pre_gz" (for *fq.gz files)

d Run mapping by typing "sh run.sh RUN_mm" (for mm10 mapping); "RUN_hs" not added yet, or just replace the ref_path in filehead of "run.sh"


 This pipelines will give cell-gene/cell-bin count matrices for each sub-library.
