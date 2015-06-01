#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use File::Basename;
use lib dirname $0;
use Sequence;
use List::Util qw (sum max);
use Bio::SeqIO;
# argumentes
my ($fasta,$type,$include);
my $aaType = $Sequence::aa;
GetOptions(
  "fasta|f:s"  =>  \$fasta,
  "type|t:s"   =>  \$type,
  "help|h"     =>  sub {&usage;exit(-1);}
);
unless ($fasta){print STDERR "Pls specify --fasta|-f\n";exit(-1);}
if(!defined($type)){
  print STDERR "Pls specify --type|-t\n"; exit(-1);
}elsif($type ne "T" && $type ne "P"){
  print STDERR "Invalid parameter for --type|-t; only i) T: Nucleotide Sequence; and ii) P: AA Sequence are legal.\n";
  exit(-1);
}
# working part
my $in = Bio::SeqIO->new(-file => $fasta, -format => 'fasta');
while(my $seqobj = $in->next_seq){
  my $id = $seqobj->id;
  my $seq = uc($seqobj->seq);
  if ($type eq 'T'){
    my @codons =  countCodonType($seq,$id);
    # id,Fragile_Codon,Total_Codon,%
    next if @codons == 0;
    print "$id\t",sum(@codons),"\t",$#codons+1,"\t",sprintf("%.3f",sum(@codons)/($#codons+1)),"\n";	
  }elsif($type eq 'P'){
    my @aa = countAAType($seq,$id);
    # id,Fragile,Facultative,Robust,%,%,%
    next if sum(@aa) == 0;
    print "$id\t",join("\t",@aa),"\t",join("\t",map{sprintf("%.3f",$_/sum(@aa))}@aa),"\n";
  }
}
# sub
sub countCodonType{
  my ($SEQ,$id) = @_;
  my @CODONTYPE;
  foreach(my $I=0;$I<length($SEQ);$I+=3){
    my $ID = substr ($SEQ,$I,3);
    # irregular CODONs and the Stop Codon was not counted;
    if ($ID !~ /^[ATCG]{3,3}$/){
    }elsif($ID =~ /TAA/ || $ID =~ /TAG/ || $ID =~ /TGA/){
      if(substr($SEQ,$I+3,3)){
        print STDERR "$id\tWarnings: There is PTC in your sequence!\n";
      }
    }else{
      push @CODONTYPE,($Sequence::codonType{$ID}=~/Fragile/?1:0);
    }
  }
  return(@CODONTYPE);
}
sub countAAType{
  my ($SEQ,$id) = @_;
  my ($AATYPE_FR,$AATYPE_FA,$AATYPE_RO) = (0,0,0);
  foreach(my $I=0;$I<length($SEQ);$I++){
    my $ID = substr($SEQ,$I,1);
    if($ID !~ /[$aaType]/ && $ID !~ /[*]/){
      print STDERR "$id\tWarnings: There are irregular AA In your sequence!\n";
    }elsif($ID =~ /[*]/){
      if (substr($SEQ,$I+1,1)){
        print STDERR "$id\tWarnings: There is Premature Stop Codon In your sequence!\n";
      }
    }else{
      $AATYPE_FR++ if $Sequence::aaType{"$ID"} eq 'Fragile';
      $AATYPE_FA++ if $Sequence::aaType{"$ID"} eq 'Facultative';
      $AATYPE_RO++ if $Sequence::aaType{"$ID"} eq 'Robust';
    }
  }
  return ($AATYPE_FR,$AATYPE_FA,$AATYPE_RO);
}
sub usage {
print STDERR <<HELP
Usage: perl $0 --fasta|-f <file.fasta> --type|-t <P/T> --help|-h
  --fasta|-f  AA / nucleotide sequences in .fasta format
  --type|-t   type of your sequence [T:nucleotide;P:protein]
  --help|-h   print this help information
Release Note:
  Function: count codon/AA type (fragile, robust or facultative); PMID: 22022272
  \@version 1.0
    Beta version.
    Irregular CODONs/AAs and the Stop Codon/* was not counted
  \@version 1.1
    Record the id with irregular CODONs/AAs
HELP
}
