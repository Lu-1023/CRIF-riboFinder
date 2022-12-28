#!/usr/bin/perl
 use strict;
 use warnings;
 #use Number::Format qw(round);
 use Getopt::Long;
 use List::MoreUtils ':all';

=head1
The scritp is used to calculate the codon motif position .the codon from frame-shifting fold change is equal or more than 35.
=cut

my $species = "";
my $frame_shift = "";
GetOptions 
(
  'sp=s'=>\$species	, # species you set 
  'fs=s'=>\$frame_shift  
 );

if ($species eq "" | $frame_shift eq "") {
	die "Missing one or more input file name(s)!
	use of script:
	perl codon_motif_pos.pl -sp <species> -fs <YES , NO , LOW,MEDIAN or ALL>) \n"} # check parameters



my %can_cod;
if ($frame_shift eq "YES"){
    #FoldChange > 35  : the high frame-shifting codons
    %can_cod = (
	    "TTC" => "", 
	    "TAT" => "",
	    "ATC" => "",
	    "TAC" => "",
        "TGG" => "",
        "AAA" => "",
        "ATA" => "",
        "ATT" => "",
        "TTT" => "",    
    );
}elsif($frame_shift eq "NO" ){
    #Foldchange = 1 and 2 :the non frame-shifting codons
    %can_cod = (
        "GCC" => "",
        "GAG" => "",
        "GAA" => "",
        "GGC" => "",
        "TCA" => "",
        "TCG" => "",
        "GCT" => "",
        "GCA"=> "",
        "ACA"=> "",
        "TCC"=> "",
        "CCA"=> "",
        "CAA"=> "",
        "ACC"=> "",
        "GGT"=> ""
    );

}elsif($frame_shift eq "MEDIAN"){
    #10 <Foldchange <=  35
    %can_cod=(
       "CAC"=>"",
       "CGT"=>"",
       "TGC" =>"",
        "GTT"=>"",
        "AGC"=>"",
        "AGA"=>"",
        "CTG"=>"",
        "CGC"=>"",
        "CTT"=>"",
        "AGG"=>"",
        "CGA"=>"",
        "CGG"=>""
    )
}elsif($frame_shift eq "LOW"){
    #2 < Foldchange <= 10 
    %can_cod=(
        "AGT"=>"",
        "ACG"=>"",
        "GGA"=>"",
        "GCG"=>"",
        "CCT"=>"",
        "TTA"=>"",
        "TCT"=>"",
        "ACT"=>"",
        "AAT"=>"",
        "GTA"=>"",
        "TTG"=>"",
        "GAC"=>"",
        "TGT"=>"",
        "AAC"=>"",
        "GGG"=>"",
        "GTG"=>"",
        "GTC"=>"",
        "CCC"=>"",
        "CTA"=>"",
        "CTC"=>"",
        "CAG"=>"",
        "AAG"=>"",
        "CCG"=>"",
        "CAT"=>"",
    )
}else{
    %can_cod = (
    'TCA' => "", 'TCC' => "", 'TCG' =>  "", 'TCT' => "", 'TTC' =>  "", 'TTT' =>  "", 'TTA' =>  "", 'TTG' =>  "",
    'TAC' =>  "", 'TAT' =>  "", 'TGC' =>  "", 'TGT' => "", 'TGG' =>  "", 'CTA' =>  "", 'CTC' =>  "", 'CTG' =>  "",
    'CTT' =>  "",'CCA' =>  "",'CCC' =>  "",'CCG' =>  "", 'CCT' =>  "", 'CAC' => "",  'CAT' => "",  'CAA' => "",  
    'CAG' => "",  'CGA' =>  "",  'CGC' =>  "",  'CGG' =>  "",  'CGT' =>  "",  'ATA' =>  "",  'ATC' =>  "",  'ATT' =>  "",
  'ATG' => "",  'ACA' =>  "",  'ACC' =>  "",  'ACG' =>  "",  'ACT' =>  "",  'AAC' =>  "",  'AAT' =>  "",  'AAA' =>  "",
    'AAG' =>  "",  'AGC' =>  "",  'AGT' =>  "",  'AGA' => "",  'AGG' => "",  'GTA' =>  "",  'GTC' =>  "",  'GTG' =>  "",
    'GTT' =>  "", 'GCA' =>  "",   'GCC' =>  "",  'GCG' =>  "",  'GCT' =>  "",  'GAC' =>"",  'GAA' => "",
    'GAG' => "",  'GGA' =>  "",  'GGC' =>  "",  'GGG' =>  "",  'GGT' =>""
    );
}


mkdir "FS";
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

if  ($frame_shift eq "YES" | $frame_shift eq "NO" | $frame_shift eq "LOW"| $frame_shift eq "MEDIAN"){

    open (hand2, ">FS/$species\_fs$frame_shift\_out.txt");
    open (hand3, ">FS/$species\_fs$frame_shift\_motif.txt"); 
    open (hand4, ">FS/$species\_fs$frame_shift\_out_same.txt");
    open (hand5, ">FS/$species\_fs$frame_shift\_motif_same.txt"); 
    open (hand6, ">FS/$species\_fs$frame_shift\_out_codon_same.txt");
    open (hand7, ">FS/$species\_fs$frame_shift\_motif_codon_same.txt"); 

    my %motif;
    foreach my $id (keys %gene)   {
        my $len = length $gene{$id};
        
        for ( my $i = 100; $i < $len - 100; $i += 3)  {
            my $seq = substr($gene{$id}, $i, 3);

            if ( exists $can_cod{$seq} )  {
                my $count = 1;
        
                Search: for (my $j = $i+3; $j < $len -100; $j += 3)  {
                    my $cod = substr($gene{$id}, $j, 3);

                    if ( exists $can_cod{$cod} )  {
                        $seq .= $cod;
                        $count ++;  } 	
                    else {
                        $motif{$seq} ++;
                        if ($count >= 3)  {
                            my $end = $i + $count*3;
                            my $pep = trans_pep($seq);
                            print hand2 "$id\t$pep\_$seq\t$i\t$end\n"; 
                            my @pro=uniq(split(//,$pep));
                            my $leng=@pro;
                            if ($leng==1){
                                print hand4 "$id\t$pep\_$seq\t$i\t$end\n"; 
                            }
                            my @codon_motif;
                            for (my $z=0 ; $z < length($seq) ; $z+=3){
                                my $cm=substr($seq,$z,3);
                                push @codon_motif,$cm; }
                                my $l=uniq(@codon_motif);
                                if ($l==1){
                                    print hand6 "$id\t$pep\_$seq\t$i\t$end\n"; 
                                }
                            if ($count >=3){
                                $i =$i + ($count-1)*3;
                        }

                            }
                        last Search;  }
                }    }  }  }

    foreach my $i (keys %motif)  {
        my $pep = trans_pep($i);
        print hand3 "$i\t$pep\t$motif{$i}\n";
        my @pro=uniq(split(//,$pep));
        my $leng=@pro;
        if ($leng==1){
            print hand5 "$i\t$pep\t$motif{$i}\n"; }

        my @codon_motif;
        for (my $z=0 ; $z < length($i) ; $z+=3){
            my $cm=substr($i,$z,3);
            push @codon_motif,$cm; }
        my $l=uniq(@codon_motif);
        if ($l==1){
                print hand7 "$i\t$pep\t$motif{$i}\n"; 
                                }
        
        }		
                            
    close hand2; close hand3;close hand4; close hand5;close hand6; close hand7;}
else{
    open (hand6, ">FS/$species\_fs$frame_shift\_out_codon_same.txt");
    open (hand7, ">FS/$species\_fs$frame_shift\_motif_codon_same.txt");
    my %motif;      
    foreach my $id (keys %gene)   {
        my $len = length $gene{$id};
        for ( my $i = 100; $i < $len - 100; $i += 3)  {
            my $seq = substr($gene{$id}, $i, 3);
            if (exists $can_cod{$seq}){
                my $count=1;
                Search: for (my $j = $i+3; $j < $len -100; $j += 3)  {
                    my $cod = substr($gene{$id}, $j, 3);
                    if ( $cod eq substr($seq,0,3) )  {
                        $seq .= $cod;
                        $count ++;  
                        } 	
                    else {
                        $motif{$seq} ++;
                        if ($count >= 3)  {
                            my $end = $i + $count*3;
                            my $pep = trans_pep($seq);
                            print hand6 "$id\t$pep\_$seq\t$i\t$end\n"; 
                            }
                        if ($count >=3){
                            $i =$i + ($count-1)*3;
                        }
                        last Search;  }
                }   }   }  }
    foreach my $i (keys %motif)  {
        my $pep = trans_pep($i);
        print hand7 "$i\t$pep\t$motif{$i}\n";
        }

    close hand6;close hand7;
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
    'GCA' => 'A',    # Alanine35 
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





