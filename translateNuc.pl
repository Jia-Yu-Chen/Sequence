#!/usr/bin/perl
use warnings;
use strict;
use Bio::SeqIO;
use Bio::Tools::CodonTable;

unless (@ARGV==2){&usage;exit(-1);}
if ($ARGV[1] ne "fa" && $ARGV[1] ne "tab"){&usage;exit(-1);}
# standard codon table
my $codonTable = Bio::Tools::CodonTable->new();
$codonTable->id(1);

my $in = Bio::SeqIO->new(-file => $ARGV[0], -format => 'fasta');
while (my $seqobj = $in->next_seq){
  if($ARGV[1] eq "tab"){
    print $seqobj->id,"\t";
    print $codonTable->translate($seqobj->seq),"\n";
  }elsif($ARGV[1] eq "fa"){
    print ">",$seqobj->id,"\n";
    print $codonTable->translate($seqobj->seq),"\n";
  }
}

sub usage{
print STDERR <<HELP
Usage:
  perl $0 <input.fa> <output format>
  output format: fa or tab
Function:
  Translate nucleotide sequence, out.fa
Author:
  Jerry Chen
Last modified:
  05/02/2015
Release notes:
  \@version 1.0 two output format
HELP
}
