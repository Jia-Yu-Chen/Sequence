#!/usr/bin/perl -w
package Sequence;

use strict;
use 5.010;
require Exporter;

1;

# 0. variables. codon type, amino acid type (PMID:22022272)
our %codonType = (
  TTT=>"Robust", TTC=>"Robust", TCT=>"Robust", TCC=>"Robust",
  TAT=>"Fragile", TAC=>"Fragile", TGT=>"Fragile", TGC=>"Fragile",
  TTA=>"Fragile", TCA=>"Fragile", TAA=>"*", TGA=>"*",
  TTG=>"Fragile",TCG=>"Fragile", TAG=>"*", TGG=>"Fragile",
  CTT=>"Robust",CTC=>"Robust", CCT=>"Robust", CCC=>"Robust",
  CAT=>"Robust",CAC=>"Robust", CGT=>"Robust", CGC=>"Robust",
  CTA=>"Robust",CTG=>"Robust", CCA=>"Robust", CCG=>"Robust",
  CAA=>"Fragile",CAG=>"Fragile", CGA=>"Fragile", CGG=>"Robust",
  ATT=>"Robust",ATC=>"Robust", ACT=>"Robust", ACC=>"Robust",
  AAT=>"Robust",AAC=>"Robust", AGT=>"Robust", AGC=>"Robust",
  ATA=>"Robust",ACA=>"Robust", AAA=>"Fragile", AGA=>"Fragile",
  ATG=>"Robust",ACG=>"Robust", AAG=>"Fragile", AGG=>"Robust",
  GTT=>"Robust",GTC=>"Robust", GCT=>"Robust", GCC=>"Robust",
  GAT=>"Robust",GAC=>"Robust", GGT=>"Robust", GGC=>"Robust",
  GTA=>"Robust",GTG=>"Robust", GCA=>"Robust", GCG=>"Robust",
  GAA=>"Fragile",GAG=>"Fragile", GGA=>"Fragile", GGG=>"Robust");
our %aaType = (
  A=>"Robust",R=>"Facultative",N=>"Robust",D=>"Robust",C=>"Fragile",
  Q=>"Fragile",E=>"Fragile",G=>"Facultative",H=>"Robust",I=>"Robust",
  L=>"Facultative",K=>"Fragile",M=>"Robust",F=>"Robust",P=>"Robust",
  S=>"Facultative",T=>"Robust",W=>"Fragile",Y=>"Fragile",V=>"Robust");
our $aa =  join("",keys %aaType);

# sub 1. reverse complement sequences
sub reverseComplement(){
  my ($seq) = @_;
  $seq = join ('', reverse (split '', $seq));
  $seq =~ tr/ATCG/TAGC/;
  $seq =~ tr/atcg/tagc/;
  $seq;
}
