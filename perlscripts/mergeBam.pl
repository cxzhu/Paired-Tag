#!/usr/bin/perl
use strict;
use warnings;

open IN, $ARGV[0] or die $!;
my $outfile = <IN>;
chomp($outfile);
my %bam_list;
while(<IN>){
	chomp;
	my @tmp = split/\s+/, $_;
	$bam_list{$tmp[0]} = $tmp[1];
}
close IN;

open OUT, "|samtools view - -b > $outfile" or die $!;

my $first = "true";
foreach my $infile (sort keys %bam_list){
	my $prefix = $bam_list{$infile};
	if($first eq "true"){
		open IN, "samtools view -h $infile|" or die $!;
		$first = "false";
	}
	else{
		open IN, "samtools view $infile|" or die $!;
	}
	while(<IN>){
		chomp;
		if(m/^\@/){
			print OUT $_."\n";
		}
		else{
			my @sp = split/\s+/, $_;
			my @sp1 = split/:/, $sp[0];
			my $output = "$sp1[0]:$sp1[1]:$sp1[2]:$sp1[3]:$sp1[4]:$sp1[5]:$sp1[6]:$prefix:$sp1[7]:$sp1[8]:$sp1[9]:$sp1[10]";
			foreach my $i (1 .. $#sp){
				$output  .= "\t$sp[$i]";
			}
			print OUT $output."\n";;
		}
	}
	close IN;
}
close OUT;



my $out_sorted = substr($outfile, 0, length($outfile)-4)."_sorted.bam";
system("samtools sort $outfile -o $out_sorted &");
