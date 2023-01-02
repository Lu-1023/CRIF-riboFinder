 # **<font color="grey"><font size=30> CRIF-riboFinder</font></font>**
<font size=5><font color="grey"><p align="right">2022.12.28</p></font></font>


- [Overview](#overview)
- [System Requirements](#system-requirements)
  - [Hardware requirements](#hardware-requirements)
  - [Software requirements ](#software-requirements-)
    - [OS Requirements ](#os-requirements-)
    - [Python Dependencies](#python-dependencies)
    - [R Dependencies](#r-dependencies)
    - [Perl Dependencies](#perl-dependencies)
  - [Installation Guide: ](#installation-guide-)
- [Demo](#demo)
- [License](#license)


# <font color="GREY">Overview</font>
***
Given that ribosome profiling could accurately indicate the position of decoding center, if, for example, a codon repeat results in +1 CRIF in a given gene, we would expect to observed stronger +1 RF signals within the CRIF box than the rest regions of that gene in the ribosome profiling data. Here, ENC (effective number of codons) means the number of codons that have at least 1 matched read at +1 frame and in the corresponding regions (baseline or +1 CRIF locus). With this process, two arrays (-1 baseline and -1 CRIF locus) are generated. The +1 CRIF locus array consists of all normalized +1 RF signals within the CRIF locus, whereas +1 baseline array consists of normalized +1 RF signals excluding 10 codons from both N- and C-terminus. Statistical significance was assessed by Student’s t-test or Welch t-test. Successful CRIF event must satisfy P value < 0.05 and foldchange > 1.5. Prediction of -1 CRIF events use the same methodology.

# <font color="GREY">System Requirements</font>
***
## <font size=5>Hardware requirements</font>
`CRIF-Finder` requires only a standard computer with enough RAM to support the in-memory operations.
## <font size=5>Software requirements </font>
### <font size=3>OS Requirements </font>
CRIF-Finder is supported for Linux. The package has been tested on the following systems:
`Linux: Ubuntu 20.04`
### <font size=3>Python Dependencies</font>
```shell
scipy
getopt
sys
re
os
numpy
argparse
random
```
### <font size=3>R Dependencies</font>
```shell
qvalue
readr
ggrepel
ggplot2
ggalt
reshape2
```
### <font size=3>Perl Dependencies</font>
```shell
strict;
warnings;
Number::Format qw(round);
Data::Dumper; 
Getopt::Long;
Statistics::Basic qw(:all nofill);
```
## <font size =5>Installation Guide: </font>
All you need to do to run the script is install python and perl.In addition, these software are used in the analysis process.The installation method is as follows(Linux: Ubuntu 20.04):
conda install
```shell
wget -c https://mirrors.bfsu.edu.cn/anaconda/miniconda/Miniconda3-latest-Linux-x86_64.sh 

chmod 777 Miniconda3-latest-Linux-x86_64.sh 

bash Miniconda3-latest-Linux-x86_64.sh
```
other softwares install
```shell
conda install -c bioconda hisat2
conda install -c bioconda bowtie
conda install -c bioconda samtools
conda install -c bioconda stringtie
conda install -c bioconda bedtools
conda install -c conda-forge python
conda install -c conda-forge perl
conda install -c r r
conda install -c bioconda blast
```
The softwares we have been used on the following version:
```shell
hisat2-align-s version 2.2.1
bowtie-align-s version 1.3.0
samtools version 1.10
stringtie 2.1.4
bedtools version 2.30.0
Python 3.8.12
Perl v5.26.2
R version 4.0.5
BLASTP
```

# <font color="GREY">Demo</font>
***
A simple example of `CRIF-Finder` is at:([https://github.com/Lu-1023/CRIF_Finder/tree/main/demo](https://github.com/Lu-1023/CRIF_Finder/tree/main/demo)
# <font color="GREY">License</font>
***
This project is covered under the [Mit License](https://github.com/Lu-1023/CRIF_Finder/blob/main/LICENSE).

