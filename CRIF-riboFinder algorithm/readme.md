
# **<font color="grey"><font size=10>CRIF-riboFinder algorithm</font></font>**
<font size=5><font color="steelblue"><p align="right">2022.11.18</p></font></font>
 <font color="steelblue">(Predict codon repeats which can significantly induce frameshifting and visualize) </font>

- [1 The prediction method for CRIF   ](#1-the-prediction-method-for-crif---)
- [2 Volcano plot and adjusted P-value ](#2-volcano-plot-and-adjusted-p-value-)
- [3 Metagene analysis to CRIF locus ](#3-metagene-analysis-to-crif-locus-)
- [4 Visualization of a line graph base on above result ](#4-visualization-of-a-line-graph-base-on-above-result-)




***
##  <font size=6>1 The prediction method for CRIF   </font>
 
Predict frameshifting events caused by codon repeat.

```shell
python frame_shift_statistic_site.py [-h] -f ${frame} -s ${organism} -d ${doc}
#-h:help information.
#-f:type = str,input 1 or 2 means +1 or -1 frame.
#-s/--sp:type=str,metavar="organism",the sample is from which species.
#-d:the location where your files are.pattern like this :/media/hp/disk4/shuang/Ribo_seq/XXX")  
```


##  <font size=6>2 Volcano plot and adjusted P-value </font>

```shell
Rscript 1.R 
```

##  <font size=6>3 Metagene analysis to CRIF locus </font>
Divide multiple bins to calculate normolized ribosome footprint signals surrounding the CRIF locus.
```shell
longest_CRIF_statistic.py [-h] -f ${frame} -s ${organism} -d ${doc}
#-f:type = str,input 1 or 2 means +1 or -1 frame
#-s/--sp:type=str,metavar="organism",the sample is from which species.
#-d:the location where your files are.pattern like this :/media/hp/disk4/shuang/Ribo_seq/XXX") 
```


##  <font size=6>4 Visualization of a line graph base on above result </font>

```shell
Rscript 2.R
```
