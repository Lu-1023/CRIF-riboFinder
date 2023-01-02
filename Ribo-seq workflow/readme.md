# Ribo-seq-workflow
 **<font color="grey"><font size=10> ANALYSIS of Ribo-Seq </font></font>**
<font size=5><font color="grey"><p align="right">2020.10.27</p></font></font>
# <font color="steelblue">STEP 1  Reference Building</font>



##   <font size=4>1 processing genome fasta</font>
remove all unwanted charaters in the chromosome header of fasta file
```shell
sed -i -e 's/ .*//' genome.fa
```
## <font size =4>2 create refFlat.txt</font>
using UCSC tools gtfToGenePred ([http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64](http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64))
```shell
~/Genomes/tools/gtfToGenePred -genePredExt genes.gtf ref.txt
```
move the 12th column to first
```shell
awk 'BEGIN {FS=OFS="\t"} {$1=$12FS$1;NF--} 1' ref.txt >ref1.txt
```
remove the last 4 columns
```shell
awk 'NF{NF-=4}1' FS="\t" OFS="\t" ref1.txt > refFlat.txt
```
## <font size=4>3 build reference</font>
creat the CDS.fa like this <kbd>100nt 5'-UTR</kbd> + <kbd>CDS</kbd> + <kbd>100nt 3'-UTR</kbd>
```shell
perl cDNA/sequence_by_refFlat.pl  -g ${name}
```
creat the longest_cDNA.fa containing the longest isoforms in this file and make the bowtie/hisat2 index for the longest_cDNA 
```shell
Rscript cDNA/longest_cDNA_index.R ${name}
```
## <font size=4>4 calculate the size of genome </font>
```shell
samtools faidx genome.fa
cut -f1,2 genome.fa.fai > sizes.genome
```
## <font size=4>5 calculate CBI  </font>
```shell
Rscript cDNA/CBI.R ${name}				
perl cDNA/UTR_length.pl -g ${name}
perl cDNA/cal_CBI_CAIavg.pl -sp ${name}	
```
## <font size =4>6 creat index for genome </font>
```shell
python ~/Genomes/extract_exons.py genes.gtf >genome.exon
python ~/Genomes/extract_splice_sites.py genes.gtf >genome.ss
hisat2-bulid -p 16 genome.fa --ss ~/Genomes/${name}/Genes/genome.ss --exon ~/Genomes/${name}/Genes/genome.exon genome
```
## <font size =4>7 creat index for rRNA </font>
```shell
hisat2-bulid -p 16 ribosomes.fa rDNA
```
[Ribo_seq_step1.sh](https://github.com/Lu-1023/Ribo-seq-workflow/blob/main/Ribo_seq_step1.sh)
***
# <font color="steelblue">STEP 2 Data Processing </font>
## <font size =4>1 trim adapter </font>
You can choose the trim tools you like or use the scrip we provided.
```shell
perl ribo_code/step1_trim3.pl -ad ${ad}
#
    ${ad} provide you adapter sequence
```
## <font size =4>2 mapping to longest CDS reference or CDS reference</font>
If you have the rRNA sequence and build the responding mapping index, firstly you should remove the reads mapping to rRNA and then the remain reads will be mapped the CDS/longest CDS reference. Otherwise , you can map the reads to the CDS/longest CDS reference directly.
```shell
1 if ${rna}=="yes"
bowtie -p 16 -v 3 --un=${line}_derRNA.fq ${Genome}/${species}/ribosome/rDNA ${sam} 2>>	de_rRNA/rRNA_stats.txt > de_rRNA/${line}_rRNAalign.aln
bowtie -p 16 -S ${Genome}/${species}/${species}_ref/bowtie/${sort} ${sam} b_align/${line}.sam 2>> b_align/mapping_report_${sort}.txt
2 if ${rna}=="no"
bowtie -p 16 -S ${Genome}/${species}/${species}_ref/bowtie/${sort} ${sam} b_align/${line}.sam 2>> b_align/mapping_report_${sort}.txt
3 
samtools view -bS -@ 16 ${line}.sam > ${line}.bam  
samtools sort -@ 16 ${line}.bam -o ${line}.sort.bam  
bedtools bamtobed -i ${line}.sort.bam > ${line}.bed
sort -k1,1 -k2,2n ${line}.bed > ${line}.sort.bed
rm ${line}.sam ${line}.bam
# 
    ${sam} represents your sample, you can merge repeats to increads  the number of the reads
```
## <font size =4>3  periodicity statistic</font>
The step will statistic the periodicity based the length of reads. You can get  the distribution for the 3 kinds of frames in different reads length.
```shell
perl ../ribo_code/step3_codon_pos_length.pl
```
## <font size =4>3 mapping to genome-wide reference</font>
And then we use the data to map to the genome-wide refrence.
```shell
1 if ${rna}=="yes"
hisat2 -p 16 -x ${Genome}/${species}/Sequence/WholeGenomeFasta/hisat2/genome -U ${line}_derRNA.fq -S align2/${line}.sam 2>> align2/mapping_report_genome_wide.txt
2 if ${rna}=="no"
hisat2 -p 16 -x ${Genome}/${species}/Sequence/WholeGenomeFasta/hisat2/genome -U ${sam} -S align2/${line}.sam 2>> align2/mapping_report_genome_wide.txt
3
samtools view -bS -@ 16 ${line}.sam > ${line}.bam  
samtools sort -@ 16 ${line}.bam -o ${line}.sort.bam
cd ballgown
mkdir ${line}
stringtie -e -B -p 16 -G ${Genome}/${species}/Genes/genes.gtf -o $line/$line.FPKM.gtf ../${line}.sort.bam
```
## <font size =4>4 convert genome.bam to genome.bed</font>
Convert the .bam files mapped to the genome-wide reference to .bed based on the periodicty statistics.
```shell
python ribo_code/step16_bamtobed.py --doc ${doc} --periodicity 70
# 
    --periodicity we choose the reads whoes the proportion of any frame is more than 70%
```
## <font size =4>5 convert to the format can be visulized</font>
To visulize the data.
```shell
sort -k1,1 -k2,2n ${line}.bed > ${line}_sort.bed
bedtools genomecov -i ${line}_sort.bed -g ${Genome}/${species}/Sequence/WholeGenomeFasta/sizes.genome -strand + -bg > positive/${line}.bedgraph   #maybe you can add -d to get the single end coverage in .bedgraph
bedtools genomecov -i ${line}_sort.bed -g ${Genome}/${species}/Sequence/WholeGenomeFasta/sizes.genome -strand - -bg > negative/${line}.bedgraph
${Genome}/tools/bedGraphToBigWig positive/${line}.bedgraph ${Genome}/${species}/Sequence/WholeGenomeFasta/sizes.genome bw/${line}.positive.bw
${Genome}/tools/bedGraphToBigWig negative/${line}.bedgraph ${Genome}/${species}/Sequence/WholeGenomeFasta/sizes.genome bw/${line}.negative.bw
```
## <font size =4>6 convert longest.bed to longest_A/P.bed</font>
For every reads we choose, we only record its signal at the A position or P position.Maybe the reads screen by periodicity or not.
```shell
1
perl ../ribo_code/step4_trimbed_v2.pl -sp ${species} -p ${position} -sort $sort_file
2
perl ../ribo_code/step4_trimbed_noscreen.pl -sp ${species} -p ${position} -sort $sort_file
#   
    -sort means you mapping to which refrence, choose "longest CDS" reference or" CDS" reference
```
## <font size =4>7 convert to the format can be visulized</font>
```shell
bedtools genomecov -split -i ${line}_A.bed -bga -strand + -g ${Genome}/${species}/${species}_ref/${sort}_cDNA.chrom > ${line}_A.bedgraph
```
## <font size =4>8 calculate the frameshifting ratio</font>
You can get the relatice frameshiting ratio and the frameshifting intensity of each genes.
```shell
perl ../ribo_code/step6_codon_pos_123.pl -sp ${species}
```
## <font size =4>9 coverage statistic</font>
Calculate and plot the RPKM and the coverage.
```shell
perl ../ribo_code/step7_cov_statv2.pl -sp ${species} -d ./ -sort ${sort}
Rscript ../ribo_code/step8_seq_status.R ${species} ${sort}
```
## <font size =4>10 reads count statistic</font>
Calculate the codon accumulate reads the median reads count for each genes.
```shell
perl ../ribo_code/step10_count_reads_by_codon_m.pl -sp ${species}
```
## <font size =4>11 the reads counts statistics of start codon and stop codon</font>
Calculate the accumulation reads at the start codon and stop codon.
```shell
perl ../ribo_code/step12_re_start_stop_count.pl -sp ${species} -s ${line}
```
[Ribo_seq_step2.sh](https://github.com/Lu-1023/Ribo-seq-workflow/blob/main/Ribo_seq_step2.sh)

