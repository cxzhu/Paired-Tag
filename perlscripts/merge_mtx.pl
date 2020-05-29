#!/usr/bin/perl
use strict;
use warnings;

open IN, $ARGV[0] or die $!;
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

