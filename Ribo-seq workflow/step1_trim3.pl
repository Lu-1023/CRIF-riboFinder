use strict;
use warnings;
use Getopt::Long;

print "This script is updated on 20200914 by suang to remove the 3' linker\n";

my $min=20;
my $adapter="";
GetOptions
(
  'ad=s'=>\$adapter  #adapter sequence
);
if ($adapter eq "") {
	die "Missing one or more input file name(s)!
	use of script:
	perl step1_trim3.pl  -ad <adapter sequence>  ) \n"} # check parameters

open (fileName, "filelist.txt") or die $!;
while (<fileName>)  {
	chomp;
	#$_ =~ s/.fastq//;
	$_ =~ s/.fq.gz//;
	#$_ =~ s/_1.fq.gz//;
	print "$_\n";
	trim3($_);   }

close fileName;

####trim 3'
sub trim3                                                      {
my ($sam)=@_; 
print "processing $sam now, trimming 3' adaptor sequences\n\n";
mkdir $sam;
my $count=0; my $id=10000000; my $idfq; my $on=0; my $si=0;
my $cca=0;my %cca=();

#open (hand1,"cat $sam.fastq |") or die $!;
#open (hand1,"zcat $sam\_1.fq.gz |") or die $!;
open (hand1,"zcat $sam.fq.gz |") or die $!;
#open (hand1,"fq/$sam.fastq") or die $!;
open (hand2,">$sam/$sam\_$min.txt");
open (hand3,">$sam\_$min.fastq");

while (<hand1>)                 {
$_ =~ s/\s+$//;
$count++;
my $mod=$count % 4;

if ($mod == 1)       {
$idfq=$_; next  }

if ($mod == 2)             {
$on=0; $si=0;
if (/^(.*$adapter)/)   { #CTGTAGGCACCATCAAT  #CACTCGGGCA #[AGCT]{4}AGATCGGAAG
$_= $1;
$_ =~ s/$adapter//;   }
else                 {
$_ =~ s/(substr($adapter,0,length($adapter)-1)|substr($adapter,0,length($adapter)-2)|substr($adapter,0,length($adapter)-3)|substr($adapter,0,length($adapter)-4)|substr($adapter,0,length($adapter)-5)).*$//;
                     }
#else                 {
#$_ =~ s/(CTGTAGGCACCATCAA|CTGTAGGCACCATCA|CTGTAGGCACCATC|CTGTAGGCACCAT|CTGTAGGCACCA|CTGTAGGCACC).*$//;
#                     }

#if (/N/)             {
#next                }

if (length($_)>$min)   {
$on=1; $id++;
$si=length($_);
print hand2 $_,"\n";
print hand3 "$idfq\n";
print hand3 $_,"\n";
print hand3 "+$id\n";}     }
elsif ($mod ==0 && $on==1) {
my $qua=substr ($_,0,$si);
print hand3 $qua,"\n";     }   }
close hand1; close hand2; close hand3; 


#### uni
open (hand1,"$sam/$sam\_$min.txt");
my %uni=(); $count=0; my %sta=(); my %siz=(); my $number=0;
while (<hand1>)            {
$_ =~ s/\s+$//; $count++;
$uni{$_}++;
$sta{substr($_,0,1)}++;
$siz{length($_)}++;        }
$number=scalar keys %uni;
close hand1;
unlink "$sam/$sam\_$min.txt";

$id=10000000;
open (hand1,">$sam/$sam\_$min\_uni.txt");
foreach my $seq (sort {$uni{$b}<=>$uni{$a}} keys %uni)  {
$id++;
print hand1 ">$id\_x$uni{$seq}\n$seq\n";                }
close hand1; %uni=(); 

open (hand1,">$sam/sum\_$sam.txt");
print hand1 "$sam\n\nTotal number of $min nt or longer reads\t$count\n";
print hand1 "Total unique species of $min nt or longer\t$number\n";
close hand1;

open (hand1,">>$sam/sum\_$sam.txt");
print hand1 "\nstart:\n";
foreach my $st (keys %sta)           {
print hand1 $st,"\t",$sta{$st},"\n"; }
%sta=();

print hand1 "\nsize:\n";
foreach my $si (sort {$a<=>$b} keys %siz)   {
print hand1 $si,"\t",$siz{$si},"\n";        }

close hand1; %siz=(); %sta=();                                   }
