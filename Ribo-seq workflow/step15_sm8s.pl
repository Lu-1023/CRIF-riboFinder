#!/usr/bin/perl
 use strict;
 use warnings;
 use Number::Format qw(round);
 use Data::Dumper; 
 use Getopt::Long;
 use Statistics::Basic qw(:all nofill);

my ($species, $flank, $coverage_percent, $codon_cov_number,$RNA_sequencing_status,$bin_number);
$species = "";
$flank = 100; # flank can not be 0
$coverage_percent = 10;
$codon_cov_number = 100;
# whether RNA seq is avaiable
 
$bin_number = 30;

GetOptions 
(
  'flank=s'=>\$flank,	# file containing peak, must have its path.
  'sp=s'=>\$species,		# name of sample

 );

if ($species eq "") {
	die "Missing one or more input file name(s)!
	use of script:
	perl step15_sm8s.pl  -sp <species name>  ) \n"} # check parameters



# This script is created by Yunkun to calculate the coverage and coverage depth of genes for ribosomal profiling. It will calculate those value for both mRNA and RPF. It will also calculate the median value of RPF, which will be used for pausing site calling.

mkdir 'sms_stat';

my %si;
open (chrom1, "/media/hp/disk1/song/Genomes/$species/$species\_ref/longest_cDNA.chrom") or die $!;

while (<chrom1>) {
	chomp;
	my @ge = split /\t/;
	$si{$ge[0]} = $ge[1];
	#print "$ge[0]\_$ge[1]\_\n";
	}
close chrom1;
	


open (fileName, "deadapter.txt") or die $!;
while (<fileName>)  {
	chomp;
	#next if /^RNA/;
	#$_ =~ s/^RPF_//;
	print "$_\n";
	coverage_cal($_);   }

close fileName;


sub coverage_cal                                                                               {

	my ($sam) = @_;


	print "============ now processing $sam =================\n\n";

 	my %rpf; my %mrna; my %rnaBED;

	open (hand1, "$sam\_A.bg") or die $!;
	

	while (<hand1>)                    {
		chomp;
  		my @a=split /\t/;
		for my $i ($a[1]..($a[2]-1)) {
			my $key = "$a[0]\t$i";
			$rpf{$key}=$a[3] if $a[3] > 0; 
			#print "$a[0]\t$i\t$a[3]\t$rpf{$key}\n";
			}
   	}
	close hand1;



	print "finish loading reference\n";

	my @can; #this will be set to calculate the meta data of distribution
	foreach my $id (sort keys %si)                       {
				

		my $inframe_pos = 0;
		my $inframe_count = 0;
		my $shift_plus_pos = 0;
		my $shift_plus_count = 0;
		my $shift_minus_pos = 0;
		my $shift_minus_count = 0;

		# count RPF in ORF
		for (my $i=$flank; $i<($si{$id}-$flank); $i+=3)              {
			# RPF inframe
			if (exists $rpf{"$id\t$i"}) {
				$inframe_pos++;                
				$inframe_count += $rpf{"$id\t$i"};
			}
			# RPF +1 shift
			my $j = $i + 1;
			if (exists $rpf{"$id\t$j"}) {
				$shift_plus_pos++;                
				$shift_plus_count += $rpf{"$id\t$j"};
			}
			# RPF -1 shift
			my $k = $i - 1;
			if (exists $rpf{"$id\t$k"}) {
				$shift_minus_pos++;                
				$shift_minus_count += $rpf{"$id\t$k"};
			}
			
		}
		my $in = round($inframe_pos/($si{$id}-200)*100, 2);
		my $plus = round($shift_plus_pos/($si{$id}-200)*100, 2);
		my $minus = round($shift_minus_pos/($si{$id}-200)*100, 2);

		push @can, $id if $in >= $coverage_percent or $inframe_pos >= $codon_cov_number; 


		next if $inframe_pos < 1;}


	meta_cov($sam, "RPF", \@can, \%rpf) ;   

}





# sub to calculate the distribution of RPF or RNA 


sub meta_cov     {

	my ($sample, $type, $candidate, $ref) = @_;
	# here, sample =  sample name, type = RPF or RNA, 
	# candidate = the gene chosen to analyze data, NOTE that selection criteria
	# ref = hash of either RPF or RNA
	

	open (hand3, ">sms_stat/$sample\_$type\_distribution_meta.txt");

	# calculate the distribution of RPF 
	# take 10 bins at UTR, ie, 20nt per bin, by default take 50 bins at gene body, normalized by ORF size

	foreach my $id ( @{$candidate} )                       {
		# first define the bin size within genes. omit those genes with ORF < 200nt.

	
		my $bin = 0;
		if (exists $si{$id} )  {		
			$bin = int(round(($si{$id}-9-200)/30,1)); print "$bin\n"; 
			 } 
		next if $bin < 6;
		print hand3 "$id";
		
		# get 5' UTR
		for my $i (0..9)  {
			my $st = $i*($flank/10);
			my $end = $st + $flank/10 - 1;
			my $tot = 0;
			for my $j ($st..$end)  {
				#print "$id\t$j\n";
				$tot += $ref->{"$id\t$j"} if exists $ref->{"$id\t$j"};
				}
			$tot=$tot/10;
			print hand3 "\t$tot";    }
    
		
		# get second codon
		my $start=0;
		my $z=$flank+3;
		$start =$ref->{"$id\t$z"} if exists $ref->{"$id\t$z"};
		
		
		print hand3 "\t$start";    
	
		# get internal 
		for my $i ( 0..($bin_number-1) )  {
			my $st = $i*$bin + $flank+6;
			my $end = $st + $bin - 1;
			my $tot = 0;
			
			for my $j ($st..$end)  {
				$tot += $ref->{"$id\t$j"} if exists $ref->{"$id\t$j"};
					}
			$tot=$tot/$bin;
			print hand3 "\t$tot";    } 
    

		# get stop codon
		my $end=0;
		my $y=$si{$id}-3-$flank;
		$end =$ref->{"$id\t$y"} if exists $ref->{"$id\t$y"};
		print hand3 "\t$end";

		# get 3' UTR
		for my $i (0..9)  {
			my $st = $i*($flank/10) + $si{$id} - $flank;
			my $end = $st + $flank/10 - 1;
			my $tot = 0;
			for my $j ($st..$end)  {
				$tot += $ref->{"$id\t$j"} if exists $ref->{"$id\t$j"};
				}
			$tot=$tot/10;
			print hand3 "\t$tot";    }
		

		# end the calculation
		print hand3 "\n";    }

	close hand3;    }					

		
	

 
