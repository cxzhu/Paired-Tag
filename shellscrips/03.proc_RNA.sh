s="sample-prefix"
mm10="path-to-mm10-STAR-ref"
mm10_rna="path-to-your-mm10-annotation-reference"
trim_galore ${s}_BC_cov.fq.gz
trim_galore -a AAAAAAAAAAAAAAAACCTGCAGGNNNNACGAATGCTCTGGCCTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN ${s}_BC_cov_trimmed.fq.gz ### trim oligo-dT primer
trim_galore -a CCTGCAGGNNNNACGAATGCTCTGGCCTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN ${s}_BC_cov_trimmed_trimmed.fq.gz ## trim N6 primer
STAR  --runThreadN 6 --genomeDir ${mm10} --readFilesIn ${s}_BC_cov_trimmed_trimmed_trimmed.fq.gz --readFilesCommand zcat --outFileNamePrefix ${s}_mm10_ --outSAMtype BAM Unsorted
samtools view -h -F 256 ${s}_mm10_Aligned.out.bam -b > ${s}\_clean.bam
samtools sort ${s}\_clean.bam ${s}_sorted
reachtools rmdup2 ${s}_sorted.bam
reachtools bam2Mtx2 ${s}_sorted_rmdup.bam ${mm10_rna}
