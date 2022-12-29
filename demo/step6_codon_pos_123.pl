#!/usr/bin/perl
 use strict;
 use warnings;
 use Number::Format qw(round);
 use Getopt::Long;

# use: creat three bedgraph files, seperating reads mapped to inframe codon, 1 position shift and 2 position shift.

# count the ratio the in-frame vs shift1 and shift2 for each gene. print out the table 

my ($species, $bin_number, $flank, $ex);
# divide each gene into 20, take 
$species = "";
$bin_number = 20;
$flank = 100;
$ex = 0; #set number of codon to exlude for analysis 

GetOptions 
(
  'b=s'=>\$bin_number,	# file containing peak, must have its path.
  'sp=s'=>\$species,		# name of sample
  'flank=s'=>\$flank,	# path and file name of genome
  'ex=s'=>\$ex,	   
 );

if ($species eq "") {
	die "Missing one or more input file name(s)!
	use of script:
	perl step4_trimbed.pl  -sp <species name>  ) \n"} # check parameters





mkdir 'frame-shift';

open (fileName, "deadapter.txt") or die $!;
while (<fileName>)  {
	chomp;
	#$_ =~ s/.fastq$//;
	next if /RNA/;
	print "$_\n";
	calculate($_);   }

close fileName;




sub calculate                                                          {

 my ($sample)=@_;
 my $ge; my %gene;  #the ncu number as key and ORF sequence as input
 my %rpf;           #coverage hash for RPF with ncu and start position as key and score as value
 my %si;            #size of each gene
 

 print "===================== now processing $sample =======================\n";

 # step1, get the length of each gene
    open (hand1, "/media/hp/disk1/song/Genomes/$species/$species\_ref/longest_cDNA.fa") or die $!;
    while (<hand1>) {
    $_=~ s/\s+$//;
    if (/^>/)       {
    $ge=$_;
    $ge=~ s/^>//;
    next;}
    $gene{$ge}.=$_; }

    close hand1;

    foreach my $j (keys %gene)  {
    $si{$j}=length($gene{$j}); 
    #print "$j\t$si{$j}\n";
                                }

 # open the reference for A site, NOTE, this file should be bedgraph and have 4 columns. 
 # generate 3 files for graphic view.

  mkdir "frame-shift/$sample"; 

  open (hand2, "$sample\_A.bg") or die $!;
  open (INFRAME, ">frame-shift/$sample/$sample\_inframe.bedgraph");
  open (s1, ">frame-shift/$sample/$sample\_shift1.bedgraph");  
  open (s2, ">frame-shift/$sample/$sample\_shift2.bedgraph");

  while (<hand2>)                    {
    chomp;
    my @a1=split /\t/;

    if ($a1[3] != 0)  {
      for my $n ($a1[1]..($a1[2]-1)) {
        $rpf{$a1[0]."\t".$n}=round($a1[3],4); } } #load hash #don't understand
  
    if ($a1[3] == 0)   {
      print INFRAME "$_\n" ; 
      print s1 "$_\n" ; 
      print s2 "$_\n" ; 
      next               }

 
    for my $p ($a1[1]..($a1[2]-1)) {
  
      my $y= ($p-$flank)%3;
      my $e=$p+1;

      if ($y==0 )   {
        print INFRAME "$a1[0]\t$p\t$e\t$a1[3]\n" ; 
        print s1 "$a1[0]\t$p\t$e\t0\n" ; 
        print s2 "$a1[0]\t$p\t$e\t0\n" ; 
               }
      elsif ($y==1)  {
        print INFRAME "$a1[0]\t$p\t$e\t0\n" ; 
        print s1 "$a1[0]\t$p\t$e\t$a1[3]\n" ;
        print s2 "$a1[0]\t$p\t$e\t0\n" ; 
               }
      elsif ($y==2)  {
        print INFRAME "$a1[0]\t$p\t$e\t0\n" ; 
        print s1 "$a1[0]\t$p\t$e\t0\n" ; 
        print s2 "$a1[0]\t$p\t$e\t$a1[3]\n" ; 
               }
    
       }  }
  close hand2; close INFRAME; close s1; close s2;
  print "finish making bg file\n";

 # count the in-frame and shift1 shift2 number
  open (hand3, ">frame-shift/$sample\_shift_count.txt") or die $!;
  print hand3 "Gene\tinframe\tshift_plus\tshift_minus\n";
  my %inframe=(); my %shift1=(); my %shift2=();

  foreach my $id (sort keys %si)                       {
    
   my $inframe=0; my $shift1=0; my $shift2=0;  
     
    for (my $i=0+$flank+$ex*3; $i<=$si{$id}-$flank-$ex*3; $i+=3)              {
       if (exists $rpf{$id."\t".$i})    {
          $inframe{$id."\t".$i}=$rpf{$id."\t".$i};
          $inframe+=$rpf{$id."\t".$i};  }             }


    for (my $i=1+$flank+$ex*3; $i<=$si{$id}-$flank-$ex*3; $i+=3)              {
       if (exists $rpf{$id."\t".$i})    {
          $shift1{$id."\t".$i}=$rpf{$id."\t".$i};
          $shift1+=$rpf{$id."\t".$i};  }             }


    for (my $i=$flank+$ex*3-1; $i<=$si{$id}-$flank-$ex*3; $i+=3)              {
       if (exists $rpf{$id."\t".$i})    {
          $shift2{$id."\t".$i}=$rpf{$id."\t".$i};
          $shift2+=$rpf{$id."\t".$i};  }             }

    print hand3 "$id\t$inframe\t$shift1\t$shift2\n";  }
  undef %rpf;
  print "finish making shift-count file\n";

 # calculate the relative position of shift
  open (INFRAME, ">frame-shift/$sample/$sample\_matrix_shift_0.txt") or die $!;
  open (s1, ">frame-shift/$sample/$sample\_matrix_shift_1.txt") or die $!;
  open (s2, ">frame-shift/$sample/$sample\_matrix_shift_2.txt") or die $!;


  foreach my $id (sort keys %si)                       {
    #print "calculate $id matrix\n";
    my $len=$si{$id};
    my @inframe = sumbin($len,$id,\%inframe);
    my @shift1 = sumbin($len,$id,\%shift1);
    my @shift2 = sumbin($len,$id,\%shift2);

    print INFRAME "$id\t";
    print s1 "$id\t";
    print s2 "$id\t";


    print INFRAME join("\t",@inframe);
    print INFRAME "\n";
    print s1 join("\t",@shift1);
    print s1 "\n";
    print s2 join("\t",@shift2);
    print s2 "\n";
                                                       }
    undef %inframe; undef %shift1; undef %shift2;
    print "finish $sample\n";

    }

sub sumbin       {
# the info: start position # region covered # $chromosome # strand # bin-size # hash for plus # hash for minus
  my ($range,$chr,$ref)=@_;
  my $count;
  my @result=();
  my $bin=int($range/$bin_number+0.5);


  for (my $i=0+$flank+$ex*3; $i<=$range-$flank-$ex*3; $i+=$bin)  {   
    $count=0;
    for (my $j=0; $j<$bin; $j++)  {
      my $pos=$j+$i; 
      $count+=$$ref{$chr."\t".$pos} if exists $$ref{$chr."\t".$pos} ;
                          }
                           
  push @result, $count;      }  
                                                             
  $ref = 0;
  return @result;    }


