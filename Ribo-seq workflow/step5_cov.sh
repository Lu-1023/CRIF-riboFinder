#!/bin/sh





### created and update on 12-31-2014 by Yunkun as a ribosomal profiling pipeline


########################################################################################
#This will be used for making the histogram of RPF and mRNA, two thing will be used separately
####################################################################################################

# Print start status message.
echo "job started"
start_time=`date +%s`


############################################################################
##make histogram for RPF
###########################################################################

species="C_elegans_Ensl_WBcel235"

cat deadapter.txt | while read line

do 
#sort -k1,1 -k2,2n $input\_A.bed > $input\_.bed

#bedtools genomecov -split -i $input/$input\_A.bed -d -strand + -g v12T0_ORF_100.chrom > $input/$input\_refA.nn.bg

#export File_name=$input"_refA"

#perl norm_bg.pl

#echo "$input for ref is done"

#this part makes the graph for IGV

bedtools genomecov -split -i ${line}_A.bed -bga -strand + -g /media/hp/disk1/song/Genomes/${species}/${species}_ref/longest_cDNA.chrom > ${line}_A.bg

done

echo "ribosomal profiling data are all done!"


##########################################################################################
###make histogram for mRNAseq
###########################################################################################


