#!/usr/bin/perl
 use strict;
 use warnings;
 use Statistics::Basic qw(:all nofill);
 use Getopt::Long;


=head1
count read density by its position on cDNA sequences. The result will be a table with every gene as row and 61 codons as columns. 
for each codon, the value will be a sum of all reads that fall into that codon. 
1, signal on initiation codon (AUG) will be count independently and named as initiator. AUG at other location will be counted independently

=cut

#set the frame to analyze 
my $frame = 0; # other option would be 1,2 or 1, -1

# set flanking sequence length of CDS
my $flank = 100;

# set region to calculate reads frequency
my $region = "M";   # chose "M" as middle to calculate the middle region, "N" to calculate front region (N terminal) and "C" to calculate rear region (C terminal)

# set width of region. If M, exclude n codon from both end. If N or C, include n codons at the given end. 
my $width = 20;
my $cutoff = $width *3*2 ;

# set species
my $species = "";

GetOptions 
(
  'flank=s'=>\$flank,	# file containing peak, must have its path.
  'sp=s'=>\$species,		# name of sample
  'region=s'=>\$region,
  'wid=s'=>\$width
 );

if ($species eq "") {
	die "Missing one or more input file name(s)!
	use of script:
	perl step10_count_reads_by_codon_m.pl  -sp <species name> ,also you can change other options ) \n"} # check parameters

mkdir "decode";

my @codon = ( 
	    'TCA' ,    # Serine
	    'TCC' ,    # Serine
	    'TCG' ,    # Serine
	    'TCT' ,    # Serine

	    'TTC' ,    # Phenylalanine
	    'TTT' ,   # Phenylalanine
	    'TTA' ,    # Leucine
	    'TTG' ,    # Leucine
	    'TAC' ,    # Tyrosine
	    'TAT' ,    # Tyrosine
	    'TAA' ,    # Stop
	    'TAG' ,    # Stop
	    'TGC',    # Cysteine
	    'TGT',    # Cysteine
	    'TGA' ,    # Stop
	    'TGG',    # Tryptophan
	    'CTA' ,    # Leucine
	    'CTC' ,    # Leucine
	    'CTG',    # Leucine
	    'CTT' ,    # Leucine

	    'CCA',    # Proline
	    'CCC' ,    # Proline
	    'CCG',    # Proline
	    'CCT' ,    # Proline

	    'CAC',    # Histidine
	    'CAT',    # Histidine
	    'CAA',    # Glutamine
	    'CAG',    # Glutamine
	    'CGA',    # Arginine
	    'CGC',    # Arginine
	    'CGG',    # Arginine
	    'CGT',    # Arginine

	    'ATA',     # Isoleucine
	    'ATC' ,    # Isoleucine
	    'ATT' ,   # Isoleucine
	    'ATG',    # Methionine

	    'ACA',    # Threonine
	    'ACC',    # Threonine
	    'ACG',    # Threonine
	    'ACT',    # Threonine

	    'AAA',    # Lysine
	    'AAG',    # Lysine
	    'AGC',    # Serine
	    'AGT',    # Serine
	    'AGA',    # Arginine
	    'AGG',    # Arginine

	    'GTA',    # Valine
	    'GTC',    # Valine
	    'GTG',    # Valine
	    'GTT',    # Valine

	    'GCA',    # Alanine
	    'GCC',    # Alanine
	    'GCG',    # Alanine
	    'GCT',    # Alanine

	    'GAC',    # Aspartic Acid
	    'GAT',    # Aspartic Acid
	    'GAA',    # Glutamic Acid
	    'GAG',    # Glutamic Acid

	    'GGA',    # Glycine
	    'GGC',    # Glycine
	    'GGG',    # Glycine
	    'GGT',    # Glycine

	    'AAT',	# Asparagine  
	    'AAC',	# Asparagine
	    );
	 

open (fileName, "deadapter.txt") or die $!;
while (<fileName>)  {
	chomp;
	
	print "$_\n";
	speed($_);   
	}

close fileName;



sub speed                                                          {

	my ($sam) = @_;

	
	print "now processing $sam\n\n";

	my $ge; my %gene;  #the ncu number as key and ORF sequence as input
	my %rpf;           #coverage hash for RPF with ncu and start position as key and score as value
		
	my %si;            #size of each gene

	open (hand1, "/media/hp/disk1/song/Genomes/$species/$species\_ref/CDS_DNA.fa") or die $!;
	while (<hand1>) {
		$_=~ s/\s+$//;

		if (/^>/)       {
    			$ge=$_;
    			$ge=~ s/^>//;
    			next;}

    			$gene{$ge}=$_; }
	close hand1;

   
 #############
#get the genome length
##############
	foreach my $j (keys %gene)  {
		$si{$j}=length($gene{$j}); 
                                }


#################
#get data from RPF bedgraph file
################# 

	open (hand2, "$sam\_A.bg") or die $!;

	while (<hand2>)                    {
		chomp;
		my @a=split /\t/;
		$rpf{$a[0]."\t".$a[1]}=$a[3] if $a[3]>0;     }
	close hand2;


###############################
#calculate coverage for each codon for a single gene
###############################

	open (hand3, ">decode/codon_reads_$sam\_$region\_$width\_frame$frame.txt");
	open (hand4, ">decode/codon_number_$sam\_$region\_$width\_frame$frame.txt");
	open (hand5, ">decode/codon_median_$sam\_$region\_$width\_frame$frame.txt");
	print hand3 "\tinitiator\t".join("\t", @codon)."\n";
	print hand4 "\t".join("\t", @codon)."\n";
	print hand5 "\tmedian\teffective_codon\ttotal_codon\n";
	foreach my $id (keys %si)   {
		
		my $init = 0; my %cod_count; my %read_count; my $e_cod=0; my $total=0; my @everycod; 

		my $m=0;
		$init = $rpf{"$id\t103"} if exists $rpf{"$id\t103"};

		next  if ($si{$id}<=$cutoff+$flank*2); 
		if ($region eq "M")  {
	
			# sum the number of codon that have reads and sum the reads for each codon
			for (my $i=$flank+$width*3+$frame; $i<($si{$id}-$flank-$width*3); $i+=3)              {
		
			if ( exists $rpf{"$id\t$i"} and $rpf{"$id\t$i"} > 0 )   {
				
				my $cod=substr($gene{$id},$i,3);
	       			$read_count{$cod} += $rpf{"$id\t$i"};
				push(@everycod,$rpf{"$id\t$i"});
				$cod_count{$cod} ++ ;
				$e_cod++;
				$total=($si{$id}-$flank*2-$cutoff)/3;			
							   }
			}
		}

		if ($region eq "C")   {
			for (my $i= $si{$id}-$flank-$width*3+$frame; $i<($si{$id}-$flank); $i+=3)              {
		
			if ( exists $rpf{"$id\t$i"} and $rpf{"$id\t$i"} > 0 )   {
				
				my $cod=substr($gene{$id},$i,3);
	       			$read_count{$cod} += $rpf{"$id\t$i"};
				push(@everycod,$rpf{"$id\t$i"}) ;
				$cod_count{$cod} ++ ;
				$e_cod ++;
				$total=($cutoff/2)/3;   }
			}
		
		}

		if ($region eq "N")   {
			for (my $i= $flank + 3 +$frame; $i< $flank + $width*3; $i+=3)              {
		
			if ( exists $rpf{"$id\t$i"} and $rpf{"$id\t$i"} > 0 )   {
				
				my $cod=substr($gene{$id},$i,3);
	       			$read_count{$cod} += $rpf{"$id\t$i"};
				push(@everycod,$rpf{"$id\t$i"}); 
				$cod_count{$cod} ++ ;
				$e_cod ++;
				$total=($cutoff/2)/3;      }
			}
		}

		$m=median(@everycod);

		# print the result 
		print hand3 "$id\t$init";
		print hand4 "$id";
		print hand5 "$id\t";
		foreach my $i (@codon) {
			if (exists $read_count{$i} )  {
				
				print hand3 "\t$read_count{$i}";
				print hand4 "\t$cod_count{$i}"; }
			else {
				print hand3 "\t0"; 
				print hand4 "\t0"; }
		}

		print hand3 "\n"; 
		print hand4 "\n";  
		print hand5 "$m\t$e_cod\t$total\n";}
	

	close hand3; close hand4;close hand5;

}
		 

