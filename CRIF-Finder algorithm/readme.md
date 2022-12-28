**<font color="grey"><font size=10>CRIF Finder </font></font>**
<font size=5><font color="steelblue"><p align="right">2022.11.18</p></font></font>
# <font color="steelblue">(Predict codon repeats which can significantly induce frameshifting and visualize) </font>



***
##  <font size=6>1 The prediction method for CRIF   </font>
 
According to the algorithm, predict significantly causes the frameshift of the codon repeat

```shell
python frame_shift_statistic_site.py [-h] -f ${frame} -s ${organism} -d ${doc}

#-f:type = str,input 1 or 2 means +1 or -1 frame
#-s/--sp:type=str,metavar="organism",the sample is from which species.
#-d:the location where your files are.pattern like this :/media/hp/disk4/shuang/Ribo_seq/XXX")  
```


##  <font size=6>2  Volcano plot and P-value adjust</font>

```shell
Rscript 1.R 
```

##  <font size=6>3 Divide multiple bins to calculate codon repeat and the situation of upstream and downstream reads   </font>

```shell
longest_CRIF_statistic.py [-h] -f ${frame} -s ${organism} -d ${doc}
#-f:type = str,input 1 or 2 means +1 or -1 frame
#-s/--sp:type=str,metavar="organism",the sample is from which species.
#-d:the location where your files are.pattern like this :/media/hp/disk4/shuang/Ribo_seq/XXX") 
```


##  <font size=6>4 Visualization of a line graph   </font>


```shell
Rscript 2.R

```
***
