#!/usr/bin/perl
use strict;
use warnings;
use Bio::Seq;
use Bio::SeqIO;

unless (@ARGV == 2){print STDERR "Usage: perl in.fa out.fa.\n";exit(-1);};

my $in  = Bio::SeqIO->new(-file => $ARGV[0], -format => 'fasta');
my $out = Bio::SeqIO->new(-format => 'Fasta',-file => ">$ARGV[1]");

while (my $seqobj = $in->next_seq){
  if( $seqobj->alphabet eq 'dna' || $seqobj->alphabet eq 'rna') {
    my $id = $seqobj->id;
    my $revcom = $seqobj->revcom;
    my $out = $out->write_seq($revcom);
  }
}
