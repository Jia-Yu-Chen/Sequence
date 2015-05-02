#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;

unless (@ARGV == 1){&usage;exit(-1);};

my $in = Bio::SeqIO->new(-file => $ARGV[0], -format => 'fasta');
while (my $seqobj = $in->next_seq) {
  my $id = $seqobj->id;
  my $seq = $seqobj->seq;
  my $length = $seqobj->length;

  if($length == 0){print $id,"\t","NA","\n";next;}

  my $gccount = @{[$seq=~m/[CGcg]/g]};
  print $id,"\t",sprintf("%.4f",$gccount/$length),"\n";
}

sub usage{
print STDERR <<HELP
Usage:
  perl $0 <input.fa>
Function:
  Report the sequence ID and GC content, tab-separated
Author:
  Jerry Chen
Last modified:
  05/02/2015
Release notes:
  \@version 1.0 beta version
HELP
}
