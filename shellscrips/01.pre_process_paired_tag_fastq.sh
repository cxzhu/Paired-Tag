s=$1 # with ${s}_R1.fq.gz and ${s}_R2.fq.gz in the working directory
p="path-to-cell-id-reference"
reachtools combine2 ${s}
zcat ${s}_combined.fq.gz | bowtie ${p} - --norc -m 1 -v 1 -S ${s}_BC.sam
reachtools convert2 ${s}_BC.sam
