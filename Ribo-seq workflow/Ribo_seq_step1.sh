#!bin/sh
###
#creat CDS_DNA.fa and longest_cDNA.fa , and hisat2/bowtie-build index longest_cDNA.fa 
#the name are same with samples in longest_cDNA_index.R and @genome in sequence_by_refFlat.pl
name="NC10"
genome="/media/hp/disk1/song/Genomes"
rna="no" #set yes or no
for i in $name
do
	cd ${genome}/${name}/Sequence/WholeGenomeFasta
	## first remove all unwanted charaters from fasta file.
	sed -i -e 's/ .*//' genome.fa
	## create own refFlat.txt with gtf
	
	cd ${genome}/${name}/Genes
	## first use gtfToGenepred
	${genome}/tools/gtfToGenePred -genePredExt genes.gtf ref.txt
	#${genome}/tools/gff3ToGenePred genes.gff ref.txt
	## next move the 12th column to first
	#awk 'BEGIN{FS=OFS="\t"} {a=$12; for (i=1;i<NF; i++) $i=$(i+1); $1=a; print}' ref.txt >ref1.txt
	awk 'BEGIN {FS=OFS="\t"} {$1=$12FS$1;NF--} 1' ref.txt >ref1.txt
	## finally remove the last 4 columns

	awk 'NF{NF-=4}1' FS="\t" OFS="\t" ref1.txt > refFlat.txt
	
	##	
	cd ${genome}
	perl cDNA/sequence_by_refFlat.pl  -g ${name}           #note: change the species in script
	Rscript cDNA/longest_cDNA_index.R ${name}             #note: change the species in script
done

#creat the genome size file
cd ${genome}/${name}/Sequence/WholeGenomeFasta
samtools faidx genome.fa
cut -f1,2 genome.fa.fai > sizes.genome

#calculate the 3'UTR and 5'UTR for the CBI
cd ${genome}

#calculate CBI_ref
Rscript cDNA/CBI.R ${name}				#note: change species in script

perl cDNA/UTR_length.pl -g ${name}

perl cDNA/cal_CBI_CAIavg.pl -sp ${name}				#note: change species in script

#creat longest index
cd ${genome}/${name}/${name}_ref
mkdir bowtie hisat2
mv *.ebwt bowtie
mv *.ht2 hisat2

#creat the genome size file and CDS_DNA index
samtools faidx CDS_DNA.fa
cut -f1,2 CDS_DNA.fa.fai > CDS_cDNA.chrom
bowtie-build --thread 16 CDS_DNA.fa CDS
mv *.ebwt bowtie

#creat genome index
cd ${genome}/${name}/Genes
python ${genome}/extract_exons.py genes.gtf >genome.exon
python ${genome}/extract_splice_sites.py genes.gtf >genome.ss
cd ${genome}/${name}/Sequence/WholeGenomeFasta
hisat2-build -p 16 genome.fa --ss ${genome}/${name}/Genes/genome.ss --exon ${genome}/${name}/Genes/genome.exon genome
mkdir hisat2
mv *.ht2 hisat2

#creat rRNA index
if [ ${rna}=="yes" ];then
	cd ${genome}/${name}/ribosome
	bowtie-build --threads 16 ribosomes.fa rDNA
fi
