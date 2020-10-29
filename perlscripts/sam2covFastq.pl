#!/usr/bin/perl
use strict;
use warnings;

open IN, "$ARGV[0]" or die $!;
my $name = substr($ARGV[0],0,length($ARGV[0])-4);
open OUT, "|gzip - > $name\_cov.fq.gz";

while(<IN>){
  next if m/^@/;
  my @tmp = split/\s+/, $_;
  next if $tmp[2] eq '*';
  my $cell_id = $tmp[2];
  my @sp = split/:/, $tmp[0];
  my $readname = "@".$sp[0].":".$cell_id.":".$sp[1];
  my $read = $sp[2];
  my $l = length($read);
  my $mark = "+";
  my $qual = substr($tmp[0], -$l, $l);
  print OUT "$readname\n$read\n$mark\n$qual\n";
}
close IN;
close OUT;
