**<font color="grey"><font size=10>Find codon repeats </font></font>**
<font size=5><font color="steelblue"><p align="right">2022.11.18</p></font></font>
<font color="steelblue">(Find all codon repeats within annotated ORF) </font>



- [1 Find the position of codon repeat motif   ](#1-find-the-position-of-codon-repeat-motif---)
- [2  Find the position of nearest stop codon at +1 or -1 frame after codon repeat](#2--find-the-position-of-nearest-stop-codon-at-1-or--1-frame-after-codon-repeat)
- [3  Merge close codon repeat ](#3--merge-close-codon-repeat-)


***
##  <font size=6>1 Find the position of codon repeat motif   </font>
 
The script is used to calculate the codon repeat position.


```shell
perl codon_motif_pos_1.pl -sp <species> -fs <YES , NO , LOW,MEDIAN or ALL>) \n"}
#-sp: type=str,the sample is from which species.
#-fs: chose fs foldchange
```


##  <font size=6>2  Find the position of nearest stop codon at +1 or -1 frame after codon repeat</font>
The frameshifting site as input and can get the distance (unit : codon) from frameshifting site to the nearest stop codon at +1 or -1 frame after codon repeat.

```shell
perl hiden_stop_codon.pl -sp <species>  -input <FS information file NAME(XXX.txt)> -output <output file NAME(XXX.txt)> ) \n"}
#-sp: type=str,the sample is from which species.
#-input:The result of the previous step
```

##  <font size=6>3  Merge close codon repeat </font>
If the distance of two codon repeats are within 10 codons(30nt),we merge two codon repeats as one.

```shell
python hsc_frame.py -sp <species> 
#-sp: type=str,the sample is from which species.

```


***
