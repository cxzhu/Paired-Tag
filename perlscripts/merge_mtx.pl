#!/usr/bin/perl
use strict;
use warnings;

### format of merge_list.txt

# prefix_of_merged_matrix
# folder_name_of_matrix_of_sub_library_1	prefix_to_add_to_cellular_barcode_of_library_1
# folder_name_of_matrix_of_sub_library_2	prefix_to_add_to_cellular_barcode_of_library_2
# ...						...

## eg
# DNA_active_merge
# CZ401_mm10_sorted_rmdup.bam_mapq10_mtx3 01
# CZ403_mm10_sorted_rmdup.bam_mapq10_mtx3 02
# CZ405_mm10_sorted_rmdup.bam_mapq10_mtx3 03
# CZ407_mm10_sorted_rmdup.bam_mapq10_mtx3 04
# CZ409_mm10_sorted_rmdup.bam_mapq10_mtx3 05
# CZ411_mm10_sorted_rmdup.bam_mapq10_mtx3 06
# CZ413_mm10_sorted_rmdup.bam_mapq10_mtx3 07
# CZ415_mm10_sorted_rmdup.bam_mapq10_mtx3 08
# CZ417_mm10_sorted_rmdup.bam_mapq10_mtx3 09
# CZ419_mm10_sorted_rmdup.bam_mapq10_mtx3 10
# CZ421_mm10_sorted_rmdup.bam_mapq10_mtx3 11
# CZ423_mm10_sorted_rmdup.bam_mapq10_mtx3 12
# CZ425_mm10_sorted_rmdup.bam_mapq10_mtx3 13
# CZ427_mm10_sorted_rmdup.bam_mapq10_mtx3 14
# CZ429_mm10_sorted_rmdup.bam_mapq10_mtx3 15
# CZ431_mm10_sorted_rmdup.bam_mapq10_mtx3 16


open IN, $ARGV[0] or die $!; ### $ARGV[0] = merge_list.txt
my %list;
my $output = <IN>;
chomp($output);
system("mkdir $output");

while(<IN>){
	chomp;
	my @tmp = split/\s+/, $_;
	$list{$tmp[0]} = $tmp[1];

}
close IN;

my %merge;
my $total = 0;
my %cell_list;
my %feat_list;
foreach my $sample (sort keys %list){
	print "Reading $sample...\n";
	my $prefix = $list{$sample};
	my %cell_idx;
	my %feat_idx;
	my $i = 0;
	open IN, "$sample/barcodes.tsv" or die $!;
	while(<IN>){
		chomp;
		$i++;
		$cell_idx{$i} = $prefix.":".$_;
	}
	close IN;
	$i = 0;
	open IN, "$sample/genes.tsv" or die $!;
	while(<IN>){
		chomp;
		$i++;
		$feat_idx{$i} = $_;
	}
	close IN;
	open IN, "$sample/matrix.mtx";
	<IN>;
	<IN>;
	<IN>;
	while(<IN>){
		chomp;
		my @tmp = split/\s+/, $_;
		my $cell_id = $cell_idx{$tmp[1]};
		my $feat_id = $feat_idx{$tmp[0]};
		my $value = $tmp[2];
		$cell_list{$cell_id} = 1;
		$feat_list{$feat_id} = 1;
		$merge{$feat_id}{$cell_id} = $value;
		$total++;
	}
	close IN;
}

### sort feat id
my $i = 0;
open OUT, ">$output/genes.tsv" or die $!;
foreach my $feat_id (sort keys %feat_list){
	$i++;
	$feat_list{$feat_id} = $i;
	print OUT $feat_id."\n";
}
close OUT;
### sort cell id
$i = 0;
open OUT, ">$output/barcodes.tsv" or die $!;
foreach my $cell_id (sort keys %cell_list){
	$i++;
	$cell_list{$cell_id} = $i;
	print OUT $cell_id."\n";
}
close OUT;

open OUT, ">$output/matrix.mtx" or die $!;
my $n_cell = keys %cell_list;
my $n_feat = keys %feat_list;
print OUT "\%\%MatrixMarket matrix coordinate real general\n\%\n$n_feat $n_cell $total\n";

foreach my $feat_id (sort keys %feat_list){
	my $fid = $feat_list{$feat_id};
	foreach my $cell_id (sort keys %{$merge{$feat_id}}){
		my $cid = $cell_list{$cell_id};
		my $value = $merge{$feat_id}{$cell_id};
		my $output = "$fid $cid $value\n";
		print OUT $output;
	}
}

close OUT;

