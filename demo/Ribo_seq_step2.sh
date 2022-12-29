##mapping , likely RNA_seq analysis
#set the species
species="hg38"
doc="/media/hp/disk1/lu/demo"
rna="yes"
position="A"
Genome="/media/hp/disk1/lu/demo"
sort="longest"   #CDS or longest
sort_file="longest_cDNA"   #or longest_cDNA 

# go to the fold storage .fq
cd ${doc}

ls *.fq.gz |sed 's/.fq.gz//'|uniq > deadapter.txt
mkdir b_align #save the result of mapping , longest_cDNA.fa as reference
mkdir de_rRNA
cat deadapter.txt | while read line 
do

	if [ ${rna} = "yes" ];then
		#remove rRNA
		echo "=============== processing $line discard rRNA ================="	
		echo "=============== processing $line =================">> de_rRNA/rRNA_stats.txt
		bowtie -p 16 -v 3 --un=${line}_derRNA.fq /media/hp/disk1/lu/demo/ribosome/rDNA ${line}.fq.gz 2>>	de_rRNA/rRNA_stats.txt > de_rRNA/${line}_rRNAalign.aln  #note ,sometimes we need "-3 4",and another we may don't need . responding to the library method 
		rm de_rRNA/${line}_rRNAalign.aln
		echo "=============== processing $line longest mRNA ================="	
		echo "=============== processing $line =================">> b_align/mapping_report_${sort}.txt
		bowtie -p 16 -S /media/hp/disk1/lu/demo/bowtie/${sort} ${line}_derRNA.fq b_align/${line}.sam 2>> b_align/mapping_report_${sort}.txt
	elif [ ${rna} = "no" ];then	
		#or  don't remove rRNA
		echo "=============== processing $line longest mRNA ================="
		echo "=============== processing $line =================">> b_align/mapping_report_${sort}.txt	
		bowtie -p 16 -S ${Genome}/${species}/${species}_ref/bowtie/${sort} ${sam} b_align/${line}.sam 2>> b_align/mapping_report_${sort}.txt
	fi
	cd b_align
	samtools view -bS -@ 16 ${line}.sam > ${line}.bam  
	echo " ===========successfully samtobam==========="
	echo " ==================sorting=================="
	samtools sort -@ 16 ${line}.bam -o ${line}.sort.bam  
	bedtools bamtobed -i ${line}.sort.bam > ${line}.bed
	sort -k1,1 -k2,2n ${line}.bed > ${line}.sort.bed
	rm ${line}.sam ${line}.bam
	cd .. 	
done
echo "       >>>>>>>>>>>>>>>>>>>>>>>2 longest finish<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"

## go to the fold storage .fq the genomewide processing need the result file from the step
cd "${doc}/b_align"
cp ../deadapter.txt ./
perl ../step3_codon_pos_length.pl
echo "      >>>>>>>>>>>>>>>>>>>>>>>3 step3_codon_pos_length finish<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"



#skip the periodicity screen

perl ../step4_trimbed_noscreen.pl -sp ${species} -p ${position} -sort $sort_file
echo "      >>>>>>>>>>>>>>>>>>>>>>>5 step4_trim bed finish<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"

cd "${doc}/b_align/noscreen"
cp ../deadapter.txt ./

cat deadapter.txt | while read line
do 
bedtools genomecov -split -i ${line}_A.bed -bga -strand + -g /media/hp/disk1/lu/demo/${sort}_cDNA.chrom > ${line}_A.bg
done
echo "      >>>>>>>>>>>>>>>>>>>>>>>5 convert finish<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"

perl ../../step6_codon_pos_123.pl -sp ${species}
echo "      >>>>>>>>>>>>>>>>>>>>>>>6 frameshift ratio finish<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"

#maybe you can also add option width to get more accurate RPKM value (skip the specific region)
perl ../../step7_cov_stat.pl -sp ${species} 
echo "      >>>>>>>>>>>>>>>>>>>>>>>7 calculate coverage finish<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"





