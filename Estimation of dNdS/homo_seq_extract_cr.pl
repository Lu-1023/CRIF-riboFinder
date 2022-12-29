use strict;
use warnings;

my %hseq; my %mseq; my $hge; my $mge;
my $file = "cr";

open (HS, "hg38.cds.fa") or die $!;


   while (<HS>)   {
       $_ =~ s/\s+$//;
       if (/^>/)   {
             $hge = $_; 
             $hge =~ s/^>//;
             next  }
        
        $hseq{$hge} .= $_;  
                     }
close HS; 


open (MM, "mus.cds.fa") or die $!;


   while (<MM>)   {
       $_ =~ s/\s+$//;
       if (/^>/)   {
             $mge = $_; 
             $mge =~ s/^>//;
             next  }
        
        $mseq{$mge} .= $_;  
	#print "$mge\n$mseq{$mge}\n";
                     }
        close MM;

 


open (hand1, "hg38_mus_$file.txt") or die $!;

open (hand2, ">hom_DNA_$file.fa");
open (hand3, ">hom_dnds_$file.csv");


while (<hand1>) {
	#print "$_";
	chomp;
	my ($qid, $sid, $piden, $qst, $qend ,$sst, $send, $qprotein, $sprotein) = split /\t/;
	#print "$qprotein\n$sprotein\n\n\n";
	my @q = split /\_/, $qid;
	my @s = split /\_/, $sid;
	$qid = "$q[0]\_$q[1]\_$q[2]";
 	$sid = "$s[0]\_$s[1]\_$s[2]";
	

	next if $piden < 80;

	#get the query sequence and add X into DNA once there are gaps
	my $qseq = "";
	$qseq = substr($hseq{$qid}, $q[5]+($qst-1)*3, ($qend-$qst+1)*3) if exists $hseq{$qid};
	#print "$qid\t$qseq\n";
	my @ins = ();
	@ins = all_match_positions("-", $qprotein);
	if (@ins)  {  #have "-"
		#print "@ins\n";
		#my $j = 0;
		foreach my $i (@ins)  {
			substr($qseq, ($i)*3, 0) = "XXX";
			 } 
	}


	#get the target sequence and add 3X into DNA once there are gaps
	my $sseq = "";
	if (exists $mseq{$sid}) {
		$sseq = substr($mseq{$sid}, $s[5]+($sst-1)*3, ($send-$sst+1)*3);}
	else {
		print "$sid sequence not exists\n"; next}

	#print "$sid\t$sseq\n";
	@ins = ();
	@ins = all_match_positions("-", $sprotein);
	if (@ins)  {  #have "-"
		#print "@ins\n";
		#my $j = 0;
		foreach my $i (@ins)  {
			#print "$i\n$sseq\n";
			substr($sseq, ($i)*3, 0) = "XXX";
			 } 
		}


	next if $qseq eq "" or $sseq eq "";
	print hand2 ">$qid\n$qseq\n";
	print hand2 ">$sid\n$sseq\n";


	my $dn = 0; my $ds = 0;
	for (my $i = 0; $i < length($qseq); $i += 3)  {
		my $qcod = substr($qseq, $i, 3);
		my $scod = substr($sseq, $i, 3);
		if ($qcod ne $scod) {

			next if $qcod eq "XXX" or $scod eq "XXX";
			my $qaa = translate($qcod);
			my $saa = translate($scod);
			
			if ($qaa eq $saa) {
				$ds += 1;}
			else {
				$dn += 1;}
		}
	}
	
	print hand3 "$qid\t$sid\t$piden\t$dn\t$ds\n";		

					}
#=cut


close hand1; close hand2; close hand3;


## this part for sub

sub all_match_positions {
    my ($regex, $string) = @_;
    my @ret;
    while ($string =~ /($regex)/g) {
        push @ret, pos($string)-length($1);
        #push @ret, [(pos($string)-length $1),pos($string)-1];  #this will report both start and end
    }
    return @ret
}


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
	#print "warning, $cod is not a codon\n";
	return "X";}
              }


