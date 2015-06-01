#!/usr/bin/perl
use warnings;
use strict;
use Bio::SeqIO;

unless (@ARGV == 2){&usage;exit(-1);}
my $in = Bio::SeqIO->new(-file => $ARGV[0], -format => 'fasta');

my $pos;
while (my $seqobj = $in->next_seq) {
  my $id  = $seqobj->id;
  my $seq = $seqobj->seq;

  if(!defined $pos){
    if ($seq =~ /$ARGV[1]/i){
	  $pos = $+[0]; # only report the first(if there are at least 2 hits
      my $count = ($seq =~ s/$ARGV[1]/$ARGV[1]/ig);
	  if ($count >= 2){
        print STDERR "WARNINGS: at least 2 hits matching $ARGV[1] in $ARGV[0].\n";
	  }
    }
    print ">$id\n",$seqobj->subseq(1,$pos),"\n";
  }else{
    print ">$id\n",substr($seq,0,$pos),"\n";
  }
}

sub usage{
print STDERR <<HELP
Usage: perl $0 <mul.fa> <up-to-and-include.seq>
Funct: 1) get the sub-seq up to <up-to-and-include.seq, eg., ATGGGA>, if more than 1 hits exist, only get the shortest sub-seq
       2) if there are multi sequences in <mul.fa>, do 1) for the first seq and trim the other seqs into sub-seq with equal length with the 1st sub-seq
HELP
}
