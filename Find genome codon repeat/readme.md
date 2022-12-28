# Find-codon-repeat
**<font color="grey"><font size=10>step1 Find codon repeats </font></font>**
<font size=5><font color="steelblue"><p align="right">2022.11.18</p></font></font>
# <font color="steelblue">(Find all codon repeats of CDS) </font>



***
##  <font size=6>1 Find the position of codon repeat motif   </font>
 
The scritp is used to calculate the codon motif position .the codon from frame-shifting fold change is equal or more than 35.


```shell
perl codon_motif_pos.pl -sp <species> -fs <YES , NO , LOW,MEDIAN or ALL>) \n"}
#-sp: type=str,the sample is from which species.
#-fs:chose fs foldchange
```


##  <font size=6>2  Count the hidden stop codon</font>
The frame-shifting site as input and can get the distance (unit : codon) from frame-shifting site to the first stop codon

```shell
perl hiden_stop_codon.pl -sp <species>  -input <FS information file NAME(XXX.txt)> -output <output file NAME(XXX.txt)> ) \n"}
#-sp: type=str,the sample is from which species.
#-input:The result of the previous step
```

##  <font size=6>3  Merge close codon repeat </font>


```shell
python hsc_frame.py -sp <species> 
#-sp: type=str,the sample is from which species.

```


***
