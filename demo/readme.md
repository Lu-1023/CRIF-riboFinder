# **<font color="grey"><font size=10> A simple example of CRIF-riboFinder </font></font>**
<font size=5><font color="grey"><p align="right">2022.12.30</p></font></font>


##   <font size=6>Overview</font>
This example shows the Ribo seq analysis of a fastq format file and the operation of the CRIF-riboFinder algorithm. Due to the limitation of file size and running time, we extracted the first 400,000 deadaptered lines of the fastq file of Ribo seq results from human liver tissue as an example (`demo.fq.gz`).
***
## <font size =6>1 Creat index for rRNA and longest</font>
<u>Note:In this script, the CDS, longest and rRNA reference used in mapping are all the demo versions extracted from the genome. All demo references(CDS_DNA.fa，longest_cDNA.fa，ribosome.fa) uploaded to GitHub([https://github.com/Lu-1023/CRIF-riboFinder/tree/main/demo](https://github.com/Lu-1023/CRIF-riboFinder/tree/main/demo))</u>:
```shell
bowtie-bulid -threads 16 ribosomes.fa rDNA
bowtie-bulid -threads 16 longest_cDNA.fa longest
```
***

##   <font size=6>2 Ribo seq workflow</font>

###   <font size=4>Run the shell script</font>

```shell
bash Ribo_seq_step2.sh
```

###   <font size=2><font color="grey">Input file</font></font>
`demo.fq.gz` as a input file

###   <font size=2><font color="grey">Expected output file</font></font>
The key results.
`b_align` is the one of output folder where we need to use the mapping results.
`de_rRNA` is the one of output folder where about the  delete rRNA results.
`b_align/noscreen` skipped the periodicity screen. In this folder,demo_A.bed records the signal at the A position.demo_A.bg is the format can be visulized.
`b_align\noscreen\coverage_stat` calculates the coverage and coverage depth of genes for ribosomal profiling. It will calculate those value for both mRNA and RPF. It will also calculate the median value of RPF, which will be used for pausing site calling.

***
## <font size =6>3 CRIF-riboFinder</font>

###  <font size=4>1 The prediction method for CRIF   </font>
 
Predict frameshifting events caused by codon repeat.

```shell
python frame_shift_statistic_site.py -f 1 -s hg38 -d /media/hp/disk1/lu/demo
python frame_shift_statistic_site.py -f 2 -s hg38 -d /media/hp/disk1/lu/demo

#-f:type = str,input 1 or 2 means +1 or -1 frame
#-s/--sp:type=str,metavar="organism",the sample is from which species.
#-d:the location where your files are.
```
###   <font size=2><font color="grey">Input file</font></font>
`CDS_DNA.fa`
`demo_A.bg`(Generated in the previous step)
`hg38_ALL_hsc_position_frame1.txt` or `hg38_ALL_hsc_position_frame2.txt`(The sites of codon repeats and the distance from frame-shifting site to the first stop codon,it's uploaded to GitHub([hg38_ALL_hsc_position_frame1.txt](https://github.com/Lu-1023/CRIF-riboFinder/tree/main/demo/b_align/noscreen) and [hg38_ALL_hsc_position_frame2.txt](https://github.com/Lu-1023/CRIF-riboFinder/tree/main/demo/b_align/noscreen))) as input files.

###   <font size=2><font color="grey">Expected output file</font></font>

`b_align\noscreen\FS` folder
##  <font size=4>2 Volcano plot and adjusted P-value</font>

```shell
Rscript 1.R s
```
###   <font size=2><font color="grey">Input file</font></font>
`b_align\noscreen\coverage_stat\demo_RPF_count.txt`(Generated in the previous step)
###   <font size=2><font color="grey">Expected output file</font></font>
The main output files
`b_align\noscreen\FS\demo\demoframe1_information_exp.csv`:The table about CRIF informations.
`b_align\noscreen\FS\demo\demo_fsgenes1_volcano.eps`:The volcano plot.


##  <font size=4>3 Metagene analysis to CRIF locus</font>
Divide multiple bins to calculate normolized ribosome footprint signals surrounding the CRIF locus.

```shell
python longest_CRIF_statistic.py -f 1 -s hg38 -d /media/hp/disk1/lu/demo
python longest_CRIF_statistic.py -f 2 -s hg38 -d /media/hp/disk1/lu/demo
#-f:type = str,input 1 or 2 means +1 or -1 frame
#-s/--sp:type=str,metavar="organism",the sample is from which species.
#-d:the location where your files are.pattern like this :/media/hp/disk4/shuang/Ribo_seq/XXX") 
```
###   <font size=2><font color="grey">Input file</font></font>
`FS/demoframe1fs_fc1.txt` or `FS/demoframe2fs_fc1.txt`(Generated in the previous step)
`hg38_ALL_hsc_position.txt`(The file containing the distance (unit : codon) from frameshifting site to the nearest stop codon at +1 or -1 frame after codon repeat. It's uploaded  to GitHub([hg38_ALL_hsc_position.txt](https://github.com/Lu-1023/CRIF-riboFinder/tree/main/demo/b_align/noscreen))

###   <font size=2><font color="grey">Expected output file</font></font>
`b_align\noscreen\FS_plot`


##  <font size=4>4 Visualization of a line graph base on above result </font>
<u>Note:In this script,because the demo data is too small to complete the last line graph. I reuploaded the original human liver FS_plot as an example of this step, and it is called FS_plot1</u>:

```shell
Rscript 2.R

```
###   <font size=2><font color="grey">Input file</font></font>
`b_align\noscreen\FS_plot1 folder`

###   <font size=2><font color="grey">Expected output file</font></font>
`liver_the reads count for repeat box1.eps`
`liver_the reads count for repeat box2.eps`
