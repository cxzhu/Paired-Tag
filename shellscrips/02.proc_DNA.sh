s="sample-prefix"
mm10="path-to-mm10-ref"
mm10_5k="path-to-your-5kb-bin-reference"
trim_galore ${s}_BC_cov.fq.gz
bowtie2 -x ${mm10} -U ${s}_BC_cov_trimmed.fq.gz --no-unal -p 8 -S ${s}_mm10.sam
samtools sort ${s}_mm10.sam ${s}_mm10_sorted.bam
reachtools rmdup2 ${s}_mm10_sorted.bam
reachtools bam2Mtx2 ${s}_mm10_sorted_rmdup.bam ${mm10_5k}
