#!/usr/bin/perl
use strict;
use warnings;


#system("mkdir tmp_list");


my $annotation = $ARGV[1]; ### annotation files with cell barcode and cell_type&state
my $bam_list = $ARGV[0]; ### bam file list contain file names and prefix


### example annotation file
### 01:03:95:01	BetaCell_WTCell_H3K27ac
### 01:03:95:04	BetaCell_KOCell_H3K9ac
###	Cell_BC			SampleInformation


### example bam list file
### SubLib1_DNA_rmdup.bam	01
### SubLib2_DNA_rmdup.bam	02
###	bam_file_name					prefix (identical to merging matrix)

## usage
## perl split_bam.pl bam_list annotation_file

open IN, $annotation or die $!;
system("mkdir split_bams");
my %cells;
my %clusters;

while(<IN>){
	chomp;
	my @tmp = split/\s+/, $_;
	my $cell_id = $tmp[0];
	my $cluster_id = $tmp[1];
	$cells{$cell_id} = $cluster_id;
	$clusters{$cluster_id} = 1;
}
close IN;

open IN, $bam_list or die $!;
my %bams;
while(<IN>){
	chomp;
	my @tmp = split/\s+/, $_;
	$bams{$tmp[0]} = $tmp[1];
}
close IN;

my %fh;

foreach my $i (sort keys %clusters){
	open $fh{$i}, "|samtools view -b - > ./split_bams/Paired-tag\_$i.bam" or die $!;
}


my $init = 0;
foreach my $bam (%bams){
	$init++;
	if($init == 1){
		open IN, "samtools view -h $bam|" or die $!;
	}
	else{
		open IN, "samtools view $bam|" or die $!;
	}
	my $prefix = $bams{$bam};

	while(<IN>){
		my $line = $_;
		if($line =~ m/^\@/){
			foreach my $i (sort keys %clusters){
				$fh{$i}->print($line);
			}
		}
		else{
			my @tmp = split/\t+/, $line;
			my $cid = substr($tmp[0], -19, 8);
			my $fcid = $prefix.":".$cid;
			next if not exists $cells{$fcid};
			my $i = $cells{$fcid};
			$fh{$i}->print($line);
		}
	}
	close IN;
}



foreach my $i (sort keys %clusters){
	close $fh{$i};
}


