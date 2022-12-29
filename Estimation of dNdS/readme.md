**<font color="grey"><font size=10>Estimation of dN/dS </font></font>**
<font size=5><font color="steelblue"><p align="right">2022.12.29</p></font></font>


# <font color="steelblue">Pipe for Estimation of dN/dS </font>

[TOC]

***
##  <font size=6>step1 Getting the pep file peptide</font>
Getting the pep file which containing the motif, type and position of codon repeat from the cds file.
```shell
perl codon_motif_pos.pl -sp ${specie} -fs ALL -group identical
perl hiden_stop_codon.pl -sp ${specie} -input ${specie}_fsALL_out_codon_same.txt -output ${specie}.txt
```


##  <font size=6>step2 Getting the Background(bg) peptide sequence</font>
The cds sequence was translated into the protein sequence as the background protein sequence.
```shell
python raw_pro_seq.py
```


##  <font size=6>3 step3 Getting the CRIF(cr) peptide sequence   </font>
Merging pep file and cds file into one file based on keyword 'id' and the 150bp downstarin of stop of codon repeat sequence or full sequence from stop of codon repeat to stop coden is obtained and translating to protein sequence according to the position of motif in pep file.

```shell
python step1_mk_mix.py
python step2_0_frame_seq.py

```

##  <font size=6>step4 BLASTP for protein sequence of bg and cr   </font>
Creating index of mus protein sequence and hg38 sequence was matched to the mus index by BLASTP.

```shell
makeblastdb -in mus_protein.fa -dbtype prot -out mus.index
blastp -db mus.index -query hg38_protein.fa -out hg38_mus.txt -outfmt '6 qseqid sseqid pident qstart qend sstart send qseq sseq' -num_threads 16 -max_hsps 1 -num_alignments 1
for i in p m
do
    for num in 150 full
    do
    makeblastdb -in mus_${i}0_protein_${num}.fa -dbtype prot -out mus_${i}_${num}.index
    blastp -db mus_${i}.index -query hg38_${i}0_protein_${num}.fa -out hg38_${specie}_${i}_${num}.txt -outfmt '6 qseqid sseqid pident qstart qend sstart send qseq sseq' -num_threads 16 -max_hsps 1 -num_alignments 1
    done
done

```

##  <font size=6>step5 Counting dnds</font>
The BLASTP result is calculated by perl script and the output(hom_dnds_bg.csv, hom_dnds_cr.csv) cotains the dn and ds value.
```shell
perl homo_seq_extract_cr.pl
```

##  <font size=6>step6 Calculating dnds and draw a box by R</font>
The final data is summarized in R and the homology which greater than 98 is screened to draw a box.
```shell
Rscript dnds_plot.R
```
