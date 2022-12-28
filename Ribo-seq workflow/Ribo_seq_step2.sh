##mapping , likely RNA_seq analysis
#set the species
species="NC12"
doc="/media/hp/disk1/song/Ribo_seq/nc"
ad="[AGCT]{4}AGATCGGAAG"
rep=2
rna="no"
position="A"
Genome="/media/hp/disk1/song/Genomes"
sort="longest"   #CDS or longest
sort_file="longest_cDNA"   #or longest_cDNA 

# go to the fold storage .fq
cd ${doc}

#deadapter
#ls *.fastq > filelist.txt
ls *.fq.gz > filelist.txt
#vi ~/Genomes/step1_trim3.pl(更改接头序列）
perl ribo_code/step1_trim3.pl -ad ${ad} ##note:change the adapter in line 40 and 42. define $min(the minist length of reads)
echo "      >>>>>>>>>>>>>>>>>>>>>>>1 deadapter finish<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"

#ls *_20.fastq |sed 's/-[12]_20.fastq//'|uniq > deadapter.txt #20 is same with $min
unset a
for i in `seq 1 ${rep}`
do
	a=${a}${i}
done
ls *_20.fastq |sed 's/-['$a']_20.fastq//'|uniq > deadapter.txt

#ls *_20.fastq |sed 's/_[12]_20.fastq//'|uniq > deadapter.txt #20 is same with $min
#ls *_20.fastq |sed 's/_20.fastq//'|uniq > deadapter.txt
mkdir b_align #save the result of mapping , longest_cDNA.fa as reference
mkdir de_rRNA
cat deadapter.txt | while read line 
do
	unset sam
	for i in `seq 1 ${rep}`
	do
		sam=${line}-${i}_20.fastq,${sam}
	done
	if [ ${rna} = "yes" ];then
		#remove rRNA
		echo "=============== processing $line discard rRNA ================="	
		echo "=============== processing $line =================">> de_rRNA/rRNA_stats.txt
		bowtie -p 16 -v 3 --un=${line}_derRNA.fq ${Genome}/${species}/ribosome/rDNA ${sam} 2>>	de_rRNA/rRNA_stats.txt > de_rRNA/${line}_rRNAalign.aln  #note ,sometimes we need "-3 4",and another we may don't need . responding to the library method 
		rm de_rRNA/${line}_rRNAalign.aln
		echo "=============== processing $line longest mRNA ================="	
		echo "=============== processing $line =================">> b_align/mapping_report_${sort}.txt
		bowtie -p 16 -S ${Genome}/${species}/${species}_ref/bowtie/${sort} ${line}_derRNA.fq b_align/${line}.sam 2>> b_align/mapping_report_${sort}.txt
	elif [ ${rna} = "no" ];then	
		# or  don't remove rRNA
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
perl ../ribo_code/step3_codon_pos_length.pl
echo "      >>>>>>>>>>>>>>>>>>>>>>>3 step3_codon_pos_length finish<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"

cd ..
mkdir align2 #save the result of mapping , genome.fa as reference
mkdir align2/ballgown
mkdir align2/positive
mkdir align2/negative
mkdir align2/bw
cat deadapter.txt | while read line
do 
	echo "=============== processing $line mapping to genome ================="	
	echo "=============== processing $line =================">> align2/mapping_report_genome_wide.txt
	if [ ${rna} = "yes" ];then
		hisat2 -p 16 -x ${Genome}/${species}/Sequence/WholeGenomeFasta/hisat2/genome -U ${line}_derRNA.fq -S align2/${line}.sam 2>> align2/mapping_report_genome_wide.txt
		# or  don't remove rRNA
	elif [ ${rna} = "no" ];then
		unset sam
		for i in `seq 1 ${rep}`
		do
			sam=${line}-${i}_20.fastq,${sam}
		done
		hisat2 -p 16 -x ${Genome}/${species}/Sequence/WholeGenomeFasta/hisat2/genome -U ${sam} -S align2/${line}.sam 2>> align2/mapping_report_genome_wide.txt
	fi
	cd align2
	samtools view -bS -@ 16 ${line}.sam > ${line}.bam  
	echo " ===========successfully samtobam==========="  
	samtools sort -@ 16 ${line}.bam -o ${line}.sort.bam
        echo " ===========successfully sorting============="
	#mapping to genomewide 
	#bedtools bamtobed -split -i ${line}.bam > ${line}_split.bed
	rm ${line}.sam 
	cd ballgown
	mkdir ${line}
	stringtie -e -B -p 16 -G ${Genome}/${species}/Genes/genes.gtf -o $line/$line.FPKM.gtf ../${line}.sort.bam
	echo " =================${line} genome finish=============="
	cd ../..
done

echo "=========== convert to bed ============="
python ribo_code/step16_bamtobed.py --doc ${doc} --periodicity 70
cat deadapter.txt | while read line
do 
	cd align2
	sort -k1,1 -k2,2n ${line}.bed > ${line}_sort.bed
	bedtools genomecov -i ${line}_sort.bed -g ${Genome}/${species}/Sequence/WholeGenomeFasta/sizes.genome -strand + -bg > positive/${line}.bedgraph   #maybe you can add -d to get the single end coverage in .bedgraph
	bedtools genomecov -i ${line}_sort.bed -g ${Genome}/${species}/Sequence/WholeGenomeFasta/sizes.genome -strand - -bg > negative/${line}.bedgraph
	${Genome}/tools/bedGraphToBigWig positive/${line}.bedgraph ${Genome}/${species}/Sequence/WholeGenomeFasta/sizes.genome bw/${line}.positive.bw
	${Genome}/tools/bedGraphToBigWig negative/${line}.bedgraph ${Genome}/${species}/Sequence/WholeGenomeFasta/sizes.genome bw/${line}.negative.bw
	cd ..
done	
	
echo "      >>>>>>>>>>>>>>>>>>>>>>>4 genome finish<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"

cd "${doc}/b_align"
perl ../ribo_code/step4_trimbed_v2.pl -sp ${species} -p ${position} -sort $sort_file
echo "      >>>>>>>>>>>>>>>>>>>>>>>5 step4_trim bed finish<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"

cat deadapter.txt | while read line
do 
bedtools genomecov -split -i ${line}_A.bed -bga -strand + -g ${Genome}/${species}/${species}_ref/${sort}_cDNA.chrom > ${line}_A.bedgraph
done
echo "      >>>>>>>>>>>>>>>>>>>>>>>5 convert finish<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"

perl ../ribo_code/step6_codon_pos_123.pl -sp ${species}
echo "      >>>>>>>>>>>>>>>>>>>>>>>6 frameshift ratio finish<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"

#maybe you can also add option width to get more accurate RPKM value (skip the specific region)
perl ../ribo_code/step7_cov_statv2.pl -sp ${species} -d ./ -sort ${sort}
echo "      >>>>>>>>>>>>>>>>>>>>>>>7 calculate coverage finish<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"

Rscript ../ribo_code/step8_seq_status.R ${species} ${sort}

perl ../ribo_code/step10_count_reads_by_codon_m.pl -sp ${species}
echo "      >>>>>>>>>>>>>>>>>>>>>>>8 count the reads on genes finish<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
cat deadapter.txt | while read line
do
        echo "${line}   ${species}" >> ./decode/detail.txt
done

#skip the periodicity screen

perl ../ribo_code/step4_trimbed_noscreen.pl -sp ${species} -p ${position} -sort $sort_file
echo "      >>>>>>>>>>>>>>>>>>>>>>>5 step4_trim bed finish<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"

cd "${doc}/b_align/noscreen"
cp ../deadapter.txt ./

cat deadapter.txt | while read line
do 
bedtools genomecov -split -i ${line}_A.bed -bga -strand + -g ${Genome}/${species}/${species}_ref/${sort}_cDNA.chrom > ${line}_A.bedgraph
done
echo "      >>>>>>>>>>>>>>>>>>>>>>>5 convert finish<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"

perl ../../ribo_code/step6_codon_pos_123.pl -sp ${species}
echo "      >>>>>>>>>>>>>>>>>>>>>>>6 frameshift ratio finish<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"

#maybe you can also add option width to get more accurate RPKM value (skip the specific region)
perl ../../ribo_code/step7_cov_statv2.pl -sp ${species} -d ./ -sort ${sort}
echo "      >>>>>>>>>>>>>>>>>>>>>>>7 calculate coverage finish<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
mkdir codon_stat
Rscript ../../ribo_code/step8_seq_status.R ${species} ${sort}

cd ${doc}/align2
cp ../deadapter.txt ./
cat deadapter.txt | while read line
do
	perl ../ribo_code/step12_re_start_stop_count.pl -sp ${species} -s ${line}
done
echo "      >>>>>>>>>>>>>>>>>>>>>>>9 count the reads of start and stop codons on genes finish<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
 
