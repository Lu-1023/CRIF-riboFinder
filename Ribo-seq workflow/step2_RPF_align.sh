#!/bin/sh





### created and update on 04-15-2014 by Yunkun as a ribosomal profiling pipeline


########################################################################################
#this is the step1, for the removal of contamination. before this, pls make sure the input file has been treated with adaptor removal
####################################################################################################

# Print start status message.
echo "job started"
start_time=`date +%s`
species="Sp_yeast_Ensl_EF2"

#defile your file name (in shell, define variables do not use $)

mkdir 'de_r+tRNA_stats'
mkdir 'bigwig'
mkdir 'bedgraph'

cat filelist.txt | while read LINE; do

	
	input=`echo "$LINE" | sed 's/.fastq$//'`
	mkdir $input
	echo "========= processing $input =================="	

	mkdir "$input/de_r+tRNA"
	# remove rRNA
	#bowtie -p 8 -v 3 --un=$input/de_r+tRNA/$input\_no_rRNA.fq ~/genomes/$species/rRNA/rDNA fq/$input\_18.fastq 2>> de_r+tRNA_stats/$input\_rRNA_stats.txt > $input/de_r+tRNA/$input\_rRNAalign.aln 

	echo "$input rRNA removal finished"

	rm $input/de_r+tRNA/$input\_rRNAalign.aln

	# remove tRNA
	#bowtie -p 8 -v 3 --un=$input/de_r+tRNA/$input\_no_rRNA_tRNA.fq ~/genomes/$species/tRNA/tRNA $input/de_r+tRNA/$input\_no_rRNA.fq 2>> de_r+tRNA_stats/$input\_tRNA_stats.txt > $input/de_r+tRNA/$input\_tRNAalign.aln 

	#echo "$input tRNA removal finished"

	rm $input/de_r+tRNA/$input\_no_rRNA.fq $input/de_r+tRNA/$input\_tRNAalign.aln 

	bowtie -p 8 -S ~/cDNA_sequence/$species\_ref/longest $input/de_r+tRNA/$input\_no_rRNA_tRNA.fq $input/$input.bt1.sam 2>> de_r+tRNA_stats/$input\_b1_stats.txt 


	# use hisat2 to do a whole genome 

	/home/yunkun/BioTools/hisat2-2.1.0/hisat2 -p 7 --dta --rna-strandness F -x /home/yunkun/genomes/$species/Sequence/hisat2/genome -U $input/de_r+tRNA/$input\_no_rRNA_tRNA.fq -S $input/HISAT_$input.sam 


	
	cd $input 

	#samtools view -Sb HISAT_$input.sam > temp.bam

	#Use samtools sort to convert the BAM file to a sorted BAM file. use 8 cores

	#samtools sort -@ 8 temp.bam -o hisat2.$input.bam

	#samtools index hisat2.$input.bam	

	rm $input.sam temp.bam


	cd ..
	### use deeptools to draw the bigwig file, for SINGLE END  #####
	# give the bw file name. Since hisat2 only give correct XS tag to strand but no take reverse complement, if use dUTP method, upper is bottom, else upper is top, --samFlagInclude 16 is taking the bottom one.since this is totally opposite, remember to exchange it if use direct ligation method.

	# Forward Strand	
	#bamCoverage -b $input/hisat2.$input.bam -o bigwig/$input.+.bw --normalizeUsing RPKM --samFlagExclude 16

	#bigWigToBedGraph bigwig/$input.+.bw bedgraph/$input.forward.bedgraph

	# Reverse Strand
	#bamCoverage -b $input/hisat2.$input.bam -o bigwig/$input.reverse.bw --normalizeUsing RPKM --samFlagInclude 16

	bigWigToBedGraph bigwig/$input.reverse.bw bedgraph/$input.reverse.bedgraph

	

	# make merged bedgraph

	cd bedgraph

	# reassign value to negative if not equal 0
	awk '{if ($4!=0){$4=-$4} print}' $input.reverse.bedgraph > $input.bot.temp.bedgraph

	sort -k1,1 -k2,2n $input.bot.temp.bedgraph > sorted.bg
	
	bedGraphToBigWig sorted.bg http://hgdownload.soe.ucsc.edu/goldenPath/$species/bigZips/$species.chrom.sizes ../bigwig/$input.-.bw
	
	cat $input.forward.bedgraph $input.bot.temp.bedgraph > $input.merge.bedgraph
	
	sort -k1,1 -k2,2n $input.merge.bedgraph > $input.sorted.bedgraph

	rm $input.forward.bedgraph $input.bot.temp.bedgraph $input.merge.bedgraph $input.reverse.bedgraph sorted.bg

	cd ../$input 


	

	samtools view -bS $input.bt1.sam > temp.bam

	#Use samtools sort to convert the BAM file to a sorted BAM file.

	samtools sort temp.bam -o $input.b1.sorted.bam

	rm $input.bt1.sam temp.bam

	# convert reads to BED format

	bedtools bamtobed -i $input.b1.sorted.bam > $input.b1.bed

	sort -k1,1 -k2,2n $input.b1.bed > $input.b1.sort.bed

	wc -l $input.b1.bed > b1.mapped.txt

	cd ..


	echo "$input alignment finished"

done


mkdir "ballgown"

# make the gtf files with estimated RNA level # 
cat filelist.txt | while read LINE; do

	input=`echo "$LINE"`
	
	echo "************** estimating RNA level of $input ********************"

	cd ballgown
 
	mkdir $input

	stringtie -e -B -p 8 -G /home/yunkun/genomes/$species/Genes/genes.gtf -o $input/$input.FPKM.gtf ../$input/hisat2.$input.bam

	cd ..
	
done


# Print end status message.
echo
echo "job finished"
end_time=`date +%s`
echo execution time was `expr $end_time - $start_time` s.
