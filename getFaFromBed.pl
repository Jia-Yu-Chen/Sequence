#!/usr/bin/perl
use warnings;
use strict;
use Bio::DB::Fasta;

my ($db,$ref,$bed);
if(@ARGV!=2){print STDERR "Usage: perl $0 ref.fa region.bed6.\n";exit(-1);}
else{open FH,$ARGV[1];}

$db = Bio::DB::Fasta->new("$ARGV[0]");
my $seq;
while(<FH>){
    chomp;
    my @t = split;
    if (defined $t[5]){
        if ($t[5] eq '+'){
            $seq = $db->seq($t[0],$t[1]+1,$t[2]);
            $seq = uc($seq);
        }else{
            $seq = $db->seq($t[0],$t[2],$t[1]+1);
            if ($t[1]+1 == $t[2]){
                $seq =~ tr/ATCGatcgn/TAGCTAGCn/;
            }
        }
    }
    print ">";
    print join("_",@t);
    print "\n$seq\n";
}
