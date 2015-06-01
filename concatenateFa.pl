#!/usr/bin/perl
use warnings;
use strict;
use Bio::SeqIO;

if (@ARGV<2){print STDERR "perl $0 1.fa 2.fa ... n.fa\n";exit(-1);}

my %outseq;

for(@ARGV){
  my $in = Bio::SeqIO->new(-file => $_, -format => 'fasta');
  while (my $seqobj = $in->next_seq) {
    my $id = $seqobj->id;
    my $seq = $seqobj->seq;
    if (exists $outseq{$id}){
      $outseq{$id} .= $seq;
    }else{
      $outseq{$id} = $seq;
    }
  }
}

foreach (keys %outseq){
  print ">$_\n",$outseq{$_},"\n";
}
