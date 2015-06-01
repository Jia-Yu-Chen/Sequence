#!/usr/bin/perl
use warnings;
use strict;
use Bio::DB::Fasta;
use Getopt::Long;

my ($db,$gpe,$nobin,$ref,$coding,$type,$intronBetweenCDSOnly);
$type = "cds";
GetOptions(
    "gpe|g:s"       =>  \$gpe,
    "nobin|b"       =>  \$nobin,
    "ref|f:s"       =>  \$ref,
    "coding|c"      =>  \$coding,
    "type|t:s"      =>  \$type,
    "cdsintron|i"   =>  \$intronBetweenCDSOnly,
    "help|h"        =>  sub{&usage;exit(-1);}
);

if(defined $gpe){open FH,$gpe;}
else{&usage();exit(-1);}

$db = Bio::DB::Fasta->new($ref);

while(<FH>){
	chomp;
	next if (/^#/);
	next if (/^chrom/);
        if ($type ne "trans" && $type ne "intron" && !defined $coding){
            print STDERR "force -c if -t !trans || !intron.\n";
            exit(-1);
        }
        $coding = 1 if ($intronBetweenCDSOnly == 1);
        
	my @field = split ("\t",$_);
	
	if (defined $nobin){unshift @field,$nobin;}
	if ($coding){
		next if /none/;
                next if $field[6] == $field[7];
	}
	my @starts = split (",",$field[9]);
	my @ends = split (",",$field[10]);
	
	my (@cds_starts,@cds_ends,@utr5_starts,@utr5_ends,@utr3_starts,@utr3_ends,@intron_starts,@intron_ends);
	my @new_starts = sort{$a<=>$b}(@starts,@field[6,7]);
	my @new_ends = sort{$a<=>$b}(@ends,@field[6,7]);
	
	if ($field[3] eq '+'){
		for(my $i=0;$i<=$#new_starts;$i++){
			if ($new_starts[$i] < $field[6]){
				push @utr5_starts,$new_starts[$i];
				push @utr5_ends,$new_ends[$i];
			}elsif($new_ends[$i] > $field[7]){
				push @utr3_starts,$new_starts[$i];
				push @utr3_ends,$new_ends[$i];
			}else{
				push @cds_starts,$new_starts[$i];
				push @cds_ends,$new_ends[$i];
			}
		}
	}
	if ($field[3] eq '-'){
		for(my $i=0;$i<=$#new_starts;$i++){
			if ($new_starts[$i] < $field[6]){
				push @utr3_starts,$new_starts[$i];
				push @utr3_ends,$new_ends[$i];
			}elsif($new_ends[$i] > $field[7]){
				push @utr5_starts,$new_starts[$i];
				push @utr5_ends,$new_ends[$i];
			}else{
				push @cds_starts,$new_starts[$i];
				push @cds_ends,$new_ends[$i];
			}
		}
	}

	my @ok;
	for(my $i=0;$i<=$#cds_starts;$i++){
		push @ok,$i if $cds_ends[$i]-$cds_starts[$i] != 0;
	}
        @cds_starts = @cds_starts[@ok];
        @cds_ends = @cds_ends[@ok];

# intron	
	if ($field[8] >= 2){
		@intron_starts = @ends[0..$field[8]-2];
		@intron_ends = @starts[1..$field[8]-1];
		if (defined $intronBetweenCDSOnly){
			@intron_starts = grep {$_ > $field[6] && $_ < $field[7]} @intron_starts;
			@intron_ends = grep {$_ > $field[6] && $_ < $field[7]} @intron_ends;
		}
	}
        
        my $outputseq;
        
        if($type eq "cds"){
            for(my $i=0; $i<=$#cds_starts; $i++){
                $outputseq .= $db->seq($field[2],$cds_starts[$i]+1,$cds_ends[$i]);
            }
        }elsif($type eq "utr5"){
            next unless defined $utr5_starts[0];
            for(my $i=0; $i<=$#utr5_starts; $i++){
                $outputseq .= $db->seq($field[2],$utr5_starts[$i]+1,$utr5_ends[$i]);
            }
        }elsif($type eq "utr3"){
            next unless defined $utr3_starts[0];
            for(my $i=0; $i<=$#utr3_starts; $i++){
                $outputseq .= $db->seq($field[2],$utr3_starts[$i]+1,$utr3_ends[$i]);
            }
        }elsif($type eq "utr"){
            if(defined $utr5_starts[0] && defined $utr3_starts[0]){
                if($utr5_starts[0] > $utr3_starts[0]){
                    for(my $i=0; $i<=$#utr5_starts; $i++){
                        $outputseq .= $db->seq($field[2],$utr5_starts[$i]+1,$utr5_ends[$i]);
                    }
                    for(my $i=0; $i<=$#utr3_starts; $i++){
                        $outputseq .= $db->seq($field[2],$utr3_starts[$i]+1,$utr3_ends[$i]);
                    }
                }else{
                    for(my $i=0; $i<=$#utr3_starts; $i++){
                        $outputseq .= $db->seq($field[2],$utr3_starts[$i]+1,$utr3_ends[$i]);
                    }
                    for(my $i=0; $i<=$#utr5_starts; $i++){
                        $outputseq .= $db->seq($field[2],$utr5_starts[$i]+1,$utr5_ends[$i]);
                    }
                }
            }elsif(defined $utr5_starts[0]){
                    for(my $i=0; $i<=$#utr5_starts; $i++){
                        $outputseq .= $db->seq($field[2],$utr5_starts[$i]+1,$utr5_ends[$i]);
                    }
            }elsif(defined $utr3_starts[0]){
                    for(my $i=0; $i<=$#utr3_starts; $i++){
                        $outputseq .= $db->seq($field[2],$utr3_starts[$i]+1,$utr3_ends[$i]);
                    }
            }else{
                next;
            }
        }elsif($type eq "intron"){
            next unless defined $intron_starts[0];
            for(my $i=0; $i<=$#intron_starts; $i++){
                $outputseq .= $db->seq($field[2],$intron_starts[$i]+1,$intron_ends[$i]);
            }
        }elsif($type eq "trans"){
            for(my $i=0; $i<=$#starts; $i++){
                $outputseq .= $db->seq($field[2],$starts[$i]+1,$ends[$i]);
            }   
        }else{
            print STDERR "Non-standard --type|-t.\n";
            exit(-1);
        }
        
        if ($field[3] eq "-"){
            $outputseq = join("",reverse(split("",$outputseq)));
            $outputseq =~ tr/atcgATCG/TAGCTAGC/;
        }else{
            $outputseq =~ tr/atcg/ATCG/;
        }
        print ">$field[12]"."|"."$field[1]\n$outputseq\n";
}

sub usage{
print STDERR <<HELP
Usage: perl $0 
    --gpe|-g:s       .gpe
    --nobin|-b       no bin column
    --ref|-f:s       ref.fa
    --coding|-c      only manipulate protein-coding genes
    --type|-t:s      type of region, force -c if -t !trans (cds,utr5,utr3,utr,intron,trans [cds])
    --cdsintron|-i   only focus on introns between coding exons when --type intron
    --help|-h        print this message

Release Note
    \@version 1.0
        get sequence from .gpe
    \@version 1.1
        add --type trans to get the sequence of the whole transcript
        force -c if -t !trans
    \@version 1.2
        utr3, utr5 or intron may be non-existance
    \@version 1.3
        force -c if -t !trans || !intron
    \@version 1.4 (01-07-2015)
        add -t utr;
        add details to examine whether the transcript is coding or not
HELP
}
