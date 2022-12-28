#!/usr/bin/perl
 use strict;
 use warnings;
 use Number::Format qw(round);

=head1
based on the ribo seq data that mapped to genome. So you'll need to 
1, mapped the RPF to genome, then use bamtobed to covert the bam to bed file 
2, based on the periodicity table( RPF mapped to longest mRNA) to trim the reads.
3, the resulted bed are pileup with bedcoverage to generate bedgraph file.

=cut

#！usr/bin/perl -w
#use List::Util qw(reduce max);
#use List::MoreUtils qw(uniq);
use strict;
use Getopt::Long;
#use List::Compare;
#use Cwd 'abs_path';

my ($sam, $species);

$sam = "";
$species = "";


GetOptions 
(
  's=s'=>\$sam,		# name of sample
  #'g=s'=>\$genome,	# path and file name of genome
  'sp=s'=>\$species,	   
 );

if ($sam eq "" or $species eq "") {
	die "Missing one or more input file name(s)!
	use of script:
	perl step12_re_start_stop_count.pl -s <sample_name> -sp <species name> ) \n"} # check parameters

#$file_in =~ s/^.+\///g;
#$refflat =~ s/^.+\///g;

# read file from bedgraph, should have 2 bedgraph for plus strand or minus strand, respectively.
my %plus; my %minus;
open (hand1, "positive/$sam.bedgraph") or die $!;
while (<hand1>)   {
	chomp;
	my ($chr, $start, $end, $val) = split /\t/;
	for my $i ($start..$end-1)    {
		$plus{$chr}{$i} = $val;   } }
close hand1;

open (hand1, "negative/$sam.bedgraph") or die $!;
while (<hand1>)   {
	chomp;
	my ($chr, $start, $end, $val) = split /\t/;
	for my $i ($start..$end-1)    {
		$minus{$chr}{$i} = $val; }  }
close hand1;

my $ge; my %genome;
open (hand2, "/media/hp/disk1/song/Genomes/$species/Sequence/WholeGenomeFasta/genome.fa") or die $!;
while (<hand2>)   {
	$_=~ s/\s+$//;
    if (/^>/)       {
    $ge=$_;
    $ge=~ s/^>//;
    next;}
    $_ =~ tr/atcg/ATCG/;
    $genome{$ge} .= $_; }

close hand2;

mkdir "start_stop";

open (REF, "/media/hp/disk1/song/Genomes/$species/Genes/refFlat.txt") or die $!;

open (Size, ">start_stop/$sam\_start_stop_counts_mRNA.txt");
print Size "gene_id\tstrand\tstart_codon\tstart_counts\tstop_codon\tstop_counts\n";

while (<REF>)   {
	#print "$_";
	chomp;
	my @a = split /\t/;
	next if !exists $genome{$a[2]}; # some genome not exists the genome
	my $start_count = 0;
	my $stop_count = 0;
	my $cod_start;
	my $cod_stop;
	my $length;
	my @e_s=split /,/,$a[9];
	my @e_e=split /,/,$a[10];
	

	

	# define CDS region, first get all exon
	my @cds;
	   
	if ($a[7] - $a[6] > 0)   {
		if ($a[3] eq "+")    { #gene in strand + 
			#Firstly, which exon is the CDS start in ?
			my $e_count=0;			
			for my $i (0..($a[8]-1)){
				if ($e_s[$i]<=$a[6] && $a[6]<=($e_e[$i]-1)){ # -1 means :end is the position+1 . For example there are 6 bases ,the start position is 0,the end position is 6 in gtf/bed/refFlat file
						$e_count=$i;
					}
				}	
			my $dvalue=$e_e[$e_count]-$a[6];
			# 3 cases, 1st : normal ; 2 nd	: Asite count likely 1st ,but codon changes ; 3 rd : Asite count and codon change
 	
			if ($dvalue >= 6){ # 1 directly take the codon from genome file and record the Asite count from the bg file
				my $p_start = $a[6]+3;
				$cod_start = substr($genome{$a[2]}, $p_start, 3);
				$start_count = $plus{$a[2]}{$p_start} if exists $plus{$a[2]}{$p_start};
				$length=$a[7]-$a[6];			
						
			}elsif ($dvalue == 5 || $dvalue==4){ # 2 codon changed : ones follow by ATG , remainings in next exon
				my $p_start = $a[6]+3;
				my $c1=$dvalue-3;
				my $c2=3-$c1;
				$cod_start = substr($genome{$a[2]}, $p_start, $c1).substr($genome{$a[2]}, $e_s[$e_count+1], $c2);
				$start_count = $plus{$a[2]}{$p_start} if exists $plus{$a[2]}{$p_start};
				$length=$a[7]-$a[6];

			}else{ # 3 codon and its count changed : codon even TG/G(start codon ATG) and codon count in next exon
				my $p_start = $e_s[$e_count+1]+3-$dvalue; #$3-d_value: get the ramains ATG (NA/TG/G) in previous exon;
				$cod_start = substr($genome{$a[2]}, $p_start, 3);
				$start_count = $plus{$a[2]}{$p_start} if exists $plus{$a[2]}{$p_start};
				$length=$a[7]-$a[6];	
			}

			#Secondly, which exon is the CDS stop in ?
			$e_count=0;			
			for my $i (0..($a[8]-1)){
				if ($e_s[$i]<=($a[7]-1) && $a[7]<=$e_e[$i]){ # -1 means :end is the position+1 . For example there are 6 bases ,the start position is 0,the end position is 6 in gtf/bed/refFlat file
						$e_count=$i;
					}
				}	
			$dvalue=$a[7]-$e_s[$e_count];
			# 2 cases, 1st : CDS stop previously only two bases(AA) or one bases(A), anther base at previous exon(U/UA) ; 2 nd : normal
 	
			if ($dvalue == 2 || $dvalue==1){ # 1 
				my $p_end = $e_e[$e_count-1]-(3-$dvalue);
				$cod_stop = substr($genome{$a[2]}, $p_end, 3-$dvalue).substr($genome{$a[2]}, $e_s[$e_count], $dvalue);
				$stop_count = $plus{$a[2]}{$p_end} if exists $plus{$a[2]}{$p_end};		
						

			}else{ # 2 directly take the codon from genome file and record the Asite count from the bg file
				my $p_end = $a[7] - 3  ;
				$cod_stop = substr($genome{$a[2]}, $p_end, 3);
				$stop_count = $plus{$a[2]}{$p_end} if exists $plus{$a[2]}{$p_end};
				
			}

		}else{ #gene in strand - 
			#Firstly, which exon is the CDS start in ?
			my $e_count=0;
			for my $i (0..($a[8]-1)){
				if ($e_s[$i]<=($a[7]-1) && $a[7]<=$e_e[$i]){ # -1 means :end is the position(or we called index ?)+1
						$e_count=$i;
					}
				}	
			my $dvalue=$a[7]-$e_s[$e_count];
			if ($dvalue >= 6){ # 1 note the - strand 
				my $p_start = $a[7] -1-3-2; #-1 means :end is the position(or we called index ?)+1 ; -3 means: remove the ATG ; -2 means : get the last position in second codon (direction : from big to small) or get the first position in second codon(direction : from small to big)
				$cod_start = substr($genome{$a[2]}, $p_start, 3);
				$cod_start = rc($cod_start);
				$start_count = $minus{$a[2]}{$p_start+2} if exists $minus{$a[2]}{$p_start+2};#correspongding line149 "-2" , codon count follow by G not the first position in second codon(direction : from small to big)
				$length=$a[7]-$a[6];			
						
			}elsif ($dvalue == 5 || $dvalue==4){ #2 codon changed : ones follow by ATG , remainings in next exon
				my $p_start = $a[7]-4;
				my $c1=$dvalue-3;
				my $c2=3-$c1;
				$cod_start = substr($genome{$a[2]}, ($e_e[$e_count-1]-$c2), $c2).substr($genome{$a[2]}, ($p_start+1-$c1), $c1);# take the codons in previous exon firstly and then take the codons follow by ATG
				$cod_start = rc($cod_start);
				$start_count = $minus{$a[2]}{$p_start} if exists $minus{$a[2]}{$p_start};
				$length=$a[7]-$a[6];
			}else{ #3 codon and its count changed : codon even GT/G(start codon GTA) and codon count in previous exon
				my $p_start = $e_e[$e_count-1]-1+$dvalue-3-2; # -1 means :end is the position(or we called index ?)+1 ; $d_value-3/-(3-$_value) : get the ramains GTA (NA/GT/G) in previous exon; -2 means : get the last position in second codon (direction : from big to small) or get the first position in second codon(direction : from small to big)
				$cod_start = substr($genome{$a[2]}, $p_start, 3);
				$cod_start = rc($cod_start);
				$start_count = $minus{$a[2]}{$p_start+2} if exists $minus{$a[2]}{$p_start+2};
				$length=$a[7]-$a[6];

			}

			#Secondly, which exon is the CDS stop in ?
			$e_count=0;			
			for my $i (0..($a[8]-1)){
				if ($e_s[$i]<=$a[6] && $a[6]<=($e_e[$i]-1)){ # -1 means :end is the position+1 . For example there are 6 bases ,the start position is 0,the end position is 6 in gtf/bed/refFlat file
						$e_count=$i;
					}
				}	
			$dvalue=$e_e[$e_count]-$a[6];
			# 2 cases, 1st : CDS stop previously only two bases(AA) or one bases(A), anther base at next exon(U/UA) ; 2 nd : normal
 	
			if ($dvalue == 2 || $dvalue==1){ # 1 
				my $p_end = $a[6];
				$cod_stop = substr($genome{$a[2]}, $p_end, $dvalue).substr($genome{$a[2]}, $e_s[$e_count+1], 3-$dvalue);
				$cod_stop = rc($cod_stop);
				$stop_count = $plus{$a[2]}{$p_end} if exists $plus{$a[2]}{$p_end};			

			}else{ # 2 directly take the codon from genome file and record the Asite count from the bg file
				my $p_end = $a[6] ;
				$cod_stop = substr($genome{$a[2]}, $a[6], 3);
				$cod_stop = rc($cod_stop);
				$stop_count = $minus{$a[2]}{$p_end} if exists $minus{$a[2]}{$p_end};
				
			}



		}   	
		
	print Size "$a[0]\_$a[1]\t$a[3]\t$cod_start\t$start_count\t$cod_stop\t$stop_count\n";	

	}
		}

close REF; close Size;	

sub rc {
	my($seq) = @_;
	my $rc = reverse $seq;
	$rc =~ tr/ATCG/TAGC/;
	return $rc ;     }

















