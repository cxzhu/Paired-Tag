argv=$1 ### argv can be QC Pre DNA RNA

### all fastq file should in name or symbol link of ${prefix}${id}_R1.fq.gz/${prefix}${id}_R2.fq.gz form
### e.g.
### CZ565_R1.fq.gz
### CZ565_R2.fq.gz
### in the $path path

path=/projects/ps-renlab/fastq/2020/2020_08_05_ChIP_A/RRPE_08052020_ChIPA/
sample_id=`seq 565 596`
sample_prefix=CZ
DNA_id="test DNA"
RNA_id="test RNA"

### set your own environment
#alias reachtools=`path to your reachtools`
mm10_bwt2_ref=/projects/ps-renlab/chz272/genome_ref/mm10/mm10
mm10_STAR_ref=/projects/ps-renlab/chz272/genome_ref/refdata-cellranger-mm10-3.0.0/star/
mm10_5k_bin_ref=/projects/ps-renlab/chz272/genome_ref/mm10.bin5k.txt
mm10_gene_matrix_ref=/projects/ps-renlab/chz272/annotations/mm10.big.txt
cell_id_ref=/projects/ps-renlab/chz272/genome_ref/PS2_CellID_ref/cell_id_full
### end of custom environment


if [[ $argv == "Pre_bz2" ]]
then
	today=`date +%Y_%m_%d`
	html_path=/projects/ps-renlab/chz272/public_html/2020/${today}/
	for q_path in $html_path 01.rawdata log
	do
		if [ ! -d "$q_path" ]
		then
			mkdir $q_path
		fi
	done
	for i in ${sample_id}
	do
		nohup sh scripts/bz2.sh $path $cell_id_ref ${sample_prefix}${i} 2>&1 > log/${sample_prefix}${i}_qc.log
	done
	mv 01.rawdata/*.html $html_path
elif [[ $argv == "Pre_gz" ]]
then
	today=`date +%Y_%m_%d`
	html_path=/projects/ps-renlab/chz272/public_html/2020/${today}/
	for q_path in $html_path 01.rawdata log
	do
		if [ ! -d "$q_path" ]
		then
			mkdir $q_path
		fi
	done
	for i in ${sample_id}
	do
		nohup sh scripts/gz.sh $path $cell_id_ref ${sample_prefix}${i} 2>&1 > log/${sample_prefix}${i}_qc.log
	done
	mv 01.rawdata/*.html $html_path
elif [[ $argv == "RUN_mm" ]]
then
	for q_path in 02.trimmed 03.mm10_mapping 04.matrices
	do
		if [ ! -d "$q_path" ]
		then
			mkdir $q_path
		fi
	done

	cd 01.rawdata
	t=0
	for i in $DNA_id
	do
		sample=${sample_prefix}${i}
		t=t+1
		if [ t>8 ]
		then
			wait
			t=0
		fi
		trim_galore ${sample}_BC_cov.fq.gz &
	done
	wait
	for i in $DNA_id
	do
		sample=${sample_prefix}${i}
		mv ${sample}_BC_cov_trimmed.fq.gz ../02.trimmed/
	done
	t=0
	for i in $RNA_id
	do
		sample=${sample_prefix}${i}
		t=t+1
		if [ t>8 ]
		then
			wait
			t=0
		fi
		trim_galore ${sample}_BC_cov.fq.gz &
	done
	wait
	t=0
	for i in $RNA_id
	do
		sample=${sample_prefix}${i}
		t=t+1
		if [ t>8 ]
		then
			wait
			t=0
		fi
		trim_galore -a AAAAAAAAAAAAAAAACCTGCAGGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN ${sample}_BC_cov_trimmed.fq.gz &
	done
	wait
	for i in $RNA_id
	do
		sample=${sample_prefix}${i}
		t=t+1
		if [ t>8 ]
		then
			wait
			t=0
		fi
		trim_galore -a CCTGCAGGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN ${sample}_BC_cov_trimmed_trimmed.fq.gz &
	done
	wait
	for i in $RNA_id
	do
		sample=${sample_prefix}${i}
		rm ${sample}_BC_cov_trimmed.fq.gz 
		mv ${sample}_BC_cov_trimmed_trimmed_trimmed.fq.gz ../02.trimmed/
	done
	for i in $DNA_id
	do
		cd ../02.trimmed/
		s=${sample_prefix}${i}
		bowtie2 -x $mm10_bwt2_ref -U ${s}_BC_cov_trimmed.fq.gz --no-unal -p 8 -S ${s}_mm10.sam
		mv ${s}_mm10.sam ../03.mm10_mapping/
		cd ../03.mm10_mapping/
		samtools sort ${s}_mm10.sam -o ${s}_mm10_sorted.bam
		reachtools rmdup2 ${s}_mm10_sorted.bam
		rm ${s}_mm10.sam
		reachtools bam2Mtx2 ${s}_mm10_sorted_rmdup.bam $mm10_5k_bin_ref
		mv ${s}*mtx2 ../04.matrices/
	done
	for i in $RNA_id
	do
		cd ../02.trimmed/
		s=${sample_prefix}${i}
		STAR  --runThreadN 8 --genomeDir $mm10_STAR_ref --readFilesIn ${s}_BC_cov_trimmed_trimmed_trimmed.fq.gz --readFilesCommand zcat --outFileNamePrefix ${s}_mm10_ --outSAMtype BAM Unsorted
		mv ${s}_mm10_Aligned.out.bam ../03.mm10_mapping/
		cd ../03.mm10_mapping
		samtools view -h -F 256 ${s}_mm10_Aligned.out.bam -b > ${s}\_clean.bam
		samtools sort ${s}\_clean.bam -o ${s}_mm10_sorted.bam
		reachtools rmdup2 ${s}_mm10_sorted.bam
		reachtools bam2Mtx2 ${s}_mm10_sorted_rmdup.bam $mm10_gene_matrix_ref
		mv ${s}*mtx2 ../04.matrices/
	done
else
	echo "== Preprocessing pipeline for Paired-seq and Paired-Tag =="
	echo "	nohup sh run.sh [ARGV] & "
	echo "		ARGV:"
	echo "		Pre_bz2		Run pre_processing from bz2 file"
	echo "		Pre_gz		Run pre_processing from gz file"
	echo "		RUN_mm		Run DNA and RNA processing pipeline, mapping to mm10"
#	echo "		RUN_hs		Run DNA and RNA processing pipeline, mapping to hg19"
fi
