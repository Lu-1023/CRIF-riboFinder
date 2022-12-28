#!/usr/bin/perl
 use strict;
 use warnings;
 #use Number::Format qw(round);
 use Getopt::Long;

=head1
you need the frame-shifting site as input and can get the distance (unit : codon) from frame-shifting site to the first stop codon

=cut
my ($species, $input, $output);

$species = "";
$input = "";
$output= "out";

mkdir "hiden_stop_codon";


GetOptions 
(
  'sp=s'=>\$species,		# 
  'input=s'=>\$input,		# the input file contains the different frame-shifting information,such as *
  'output=s'=>\$output,		# the distance (unit : codon) from frame-shifting site to the first stop codon
  
 );

if ($species eq "" or $input eq "" or $output eq "" ) {
	die "Missing one or more input file name(s)!
	use of script:
	perl hiden_stop_codon.pl -sp <species>  -input <FS information file NAME(XXX.txt)> -output <output file NAME(XXX.txt)> ) \n"} # check parameters

# load reference sequences
my %gene; my %cds; my $ge;
open (hand1, "/media/hp/disk1/song/Genomes/$species/$species\_ref/CDS_DNA.fa") or die $!;
while (<hand1>) {
    $_=~ s/\s+$//;
    if (/^>/)       {
      $ge=$_;
      $ge=~ s/^>//;
      next;}
    $gene{$ge} .= $_;
}
close hand1;

open (hand1, "FS/$input") or die $!;
open (hand2,">hiden_stop_codon/$output") ;# id motif start stop +1distance -1 distance
print hand2 "id\tmotif\tstart\tstop\t +1_stop\t-1_stop\n"; #codon pos
open (hand3,">hiden_stop_codon/pep_$output") ;# id motif start stop +1distance -1 distance
print hand3 "id\tmotif\tstart\tstop\t +1_stop\t-1_stop\t+1_pep\t-1_pep\n"; #shift pep

while (<hand1>)     {
	chomp;
	my ($id, $motif, $start, $stop) = split /\t/;
	if (exists $gene{$id}  )  {
		
		my ($shift_plus_dna, $shift_plus_pep) = shift_pep($gene{$id}, $start, 1);      # assume that shifting occurs at 3th codon of that sequence stretch
   		my ($shift_minus_dna, $shift_minus_pep) = shift_pep($gene{$id}, $stop, -1);
        my $positive_d=$start+length($shift_plus_dna)+1;
        my $negative_d=$start+length($shift_minus_dna)-1;
        $motif =~ m/(.*)_.*/;
        my $real_l=length($1)-2;
        my $up=substr($gene{$id},100,$start-100);
        my ($up_dna,$up_pep)=trans_pep($up);
        print hand2 "$id\t$motif\t$start\t$stop\t$positive_d\t$negative_d\n";
        print hand3 "$id\t$motif\t$start\t$stop\t$positive_d\t$negative_d\t$up_pep$shift_plus_pep\t$up_pep$shift_minus_pep\n";
		
	     }    }

close hand1; close hand2; 

sub shift_pep     {

    my ($seq, $start, $shift) = @_;

    my $downstream = substr($seq, $start  + $shift ); 
    my ($downdna, $downpep) = trans_pep($downstream);
    return ($downdna, $downpep);
  }


sub trans_pep     {

    my $seq = shift;

          my $pep; my $dna;
          TRANSLATE: for (my $j=0; $j < length($seq); $j+=3)     {
              my $cod = substr($seq, $j, 3);
              $dna .= $cod;
              my $aa = translate($cod);
              $pep .= $aa ;   
              last TRANSLATE if $aa eq '*';                   }
              
    return ($dna, $pep);       }




sub translate {
 my ($cod) = shift;
 
 my (%codon2aa) = (
    
    'TCA' => 'S',    # Serine
    'TCC' => 'S',    # Serine
    'TCG' => 'S',    # Serine
    'TCT' => 'S',    # Serine
    'TTC' => 'F',    # Phenylalanine
    'TTT' => 'F',    # Phenylalanine
    'TTA' => 'L',    # Leucine
    'TTG' => 'L',    # Leucine
    'TAC' => 'Y',    # Tyrosine
    'TAT' => 'Y',    # Tyrosine
    'TAA' => '*',    # Stop
    'TAG' => '*',    # Stop
    'TGC' => 'C',    # Cysteine
    'TGT' => 'C',    # Cysteine
    'TGA' => '*',    # Stop
    'TGG' => 'W',    # Tryptophan
    'CTA' => 'L',    # Leucine
    'CTC' => 'L',    # Leucine
    'CTG' => 'L',    # Leucine
    'CTT' => 'L',    # Leucine
    'CCA' => 'P',    # Proline
    'CCC' => 'P',    # Proline
    'CCG' => 'P',    # Proline
    'CCT' => 'P',    # Proline
    'CAC' => 'H',    # Histidine
    'CAT' => 'H',    # Histidine
    'CAA' => 'Q',    # Glutamine
    'CAG' => 'Q',    # Glutamine
    'CGA' => 'R',    # Arginine
    'CGC' => 'R',    # Arginine
    'CGG' => 'R',    # Arginine
    'CGT' => 'R',    # Arginine
    'ATA' => 'I',    # Isoleucine
    'ATC' => 'I',    # Isoleucine
    'ATT' => 'I',    # Isoleucine
    'ATG' => 'M',    # Methionine
    'ACA' => 'T',    # Threonine
    'ACC' => 'T',    # Threonine
    'ACG' => 'T',    # Threonine
    'ACT' => 'T',    # Threonine
    'AAC' => 'N',    # Asparagine
    'AAT' => 'N',    # Asparagine
    'AAA' => 'K',    # Lysine
    'AAG' => 'K',    # Lysine
    'AGC' => 'S',    # Serine
    'AGT' => 'S',    # Serine
    'AGA' => 'R',    # Arginine
    'AGG' => 'R',    # Arginine
    'GTA' => 'V',    # Valine
    'GTC' => 'V',    # Valine
    'GTG' => 'V',    # Valine
    'GTT' => 'V',    # Valine
    'GCA' => 'A',    # Alanine
    'GCC' => 'A',    # Alanine
    'GCG' => 'A',    # Alanine
    'GCT' => 'A',    # Alanine
    'GAC' => 'D',    # Aspartic Acid
    'GAT' => 'D',    # Aspartic Acid
    'GAA' => 'E',    # Glutamic Acid
    'GAG' => 'E',    # Glutamic Acid
    'GGA' => 'G',    # Glycine
    'GGC' => 'G',    # Glycine
    'GGG' => 'G',    # Glycine
    'GGT' => 'G',    # Glycine
    );
 
 if (exists $codon2aa{$cod}) {
  	return $codon2aa{$cod};  } 
 else { 
	print "warning, $cod is not a codon\n";
	return "X";}
              }
