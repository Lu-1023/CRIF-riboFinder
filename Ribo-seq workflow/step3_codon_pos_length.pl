use strict;
use warnings;
#use Statistics::Histogram;
use Number::Format qw/round/;

print "This script is updated on 20200916 by shuang to assess the right 5'end flank of RPF protected fragment, which determine the codon on A site\n";

my $up = 40;
my $low = 25;

mkdir 'codon_stat';

open (fileName, "deadapter.txt") or die $!;
while (<fileName>)  {
	chomp;
	print "$_\n";
	sta($_);   }

close fileName;

sub sta                                                      {

	my ($sam)=@_; 

	my %count=(); my $count=0; my %p0=(); my %pplus=(); my %pminus=(); my %size=(); my %strand=(); my %tot;

	open (hand1, "$sam.sort.bed") or die $!;

	while (<hand1>)                                             {

	$_=~ s/\s+$//;
	my @a1=split /\t/; 

	$strand{$a1[5]}++;

	if ($a1[5] eq '-') {next;}


	my $size=$a1[2]-$a1[1];

	$size{$size}++;

	if ($size<$low or $size>$up) {next;}

	$count{$size}++;
	$count++;

	         #assume A site at 16nt from 5' end
	my $f12=($a1[1]+15-100)%3;
	
	if ($f12 == 0) { 
		$tot{'0'} ++;
		$p0{$size} ++; }
	elsif ($f12 == 1) {
		$tot{'pplus'} ++;		
		$pplus{$size} ++; }
	else {
		$tot{'pminus'} ++;		
		$pminus{$size} ++;}

	                                   }
	close hand1;

	open (hand2, ">codon_stat/$sam.txt");
	open (hand3, ">>codon_stat/total.txt");

	foreach (keys %strand)  {
		print hand2 "$_\t$strand{$_}\n"; }

	foreach (sort keys %size)  {
		print hand2 "$_\t$size{$_}\n"; }

	print hand2 "The mapped reads of $low-$up nt are $count\n";


	print hand2 "length\tnumber_of_reads\t0\tplus\tminus\n";

	for my $si ($low..$up)  {

		if (exists $count{$si} and $count{$si} > 0)   {

			my $ra12=round($p0{$si}/$count{$si}*100);

			my $ra13=round($pplus{$si}/$count{$si}*100);

			my $ra14=round($pminus{$si}/$count{$si}*100);     

			print hand2 "$si\t$size{$si}\t$ra12\t$ra13\t$ra14\n" if exists $size{$si};
		
		}
                         }


	print hand3 "sample\t0\tplus\tminus\n";
	print hand3 "$sam\t$tot{'0'}\t$tot{'pplus'}\t$tot{'pminus'}\n";
	print "$sam is done\n";

	close hand2;
	close hand3;                                                          }
          


