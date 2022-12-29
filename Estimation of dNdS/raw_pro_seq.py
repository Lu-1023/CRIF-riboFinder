import pandas as pd
import re
species = ['hg38']
path = '/media/hp/disk4/YT/protein_homology/ref/CRIF/'
codonTable = {
        'AUA': 'I', 'AUC': 'I', 'AUU': 'I', 'AUG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACU': 'T',
        'AAC': 'N', 'AAU': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGU': 'S', 'AGA': 'R', 'AGG': 'R',
        'CUA': 'L', 'CUC': 'L', 'CUG': 'L', 'CUU': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCU': 'P',
        'CAC': 'H', 'CAU': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGU': 'R',
        'GUA': 'V', 'GUC': 'V', 'GUG': 'V', 'GUU': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCU': 'A',
        'GAC': 'D', 'GAU': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGU': 'G',
        'UCA': 'S', 'UCC': 'S', 'UCG': 'S', 'UCU': 'S',
        'UUC': 'F', 'UUU': 'F', 'UUA': 'L', 'UUG': 'L',
        'UAC': 'Y', 'UAU': 'Y', 'UAA': '*', 'UAG': '*',
        'UGC': 'C', 'UGU': 'C', 'UGA': '*', 'UGG': 'W'}

for specie in species:
    cds_fa = open(path + specie + '.cds.fa', 'r')
    cds_txt = open(path + specie + '.cds.txt', 'w')
    for line in cds_fa:
        if '>' in line:
            line = line.replace('>', '').replace('\n', '\t')
            cds_txt.write(line)
        else:
                cds_txt.write(line)
    cds = pd.read_csv(path + specie + '.cds.txt', delimiter='\t', header=None, names=['ID', 'SEQ'], encoding='utf-8')
    pro_f = open(path + specie + '_protein.fa', 'w')
    cds = cds.drop_duplicates()
    ID = cds['ID']
    SEQ = cds['SEQ']
    ID_v = list(map(str, ID))
    SEQ_v = list(map(str, SEQ))
    for i in range(0, len(ID_v)):
        rna = re.sub('T', 'U', SEQ_v[i][100:-100])
        proteinSeq = ""
        for codonStart in range(0, len(SEQ_v[i][100:-100]), 3):
            codon = rna[codonStart:codonStart + 3]
            if codon in codonTable:
                proteinSeq += codonTable[codon]
        #if len(proteinSeq) >= 100:
        pro_f.write('>' + ID_v[i] + '\n' + proteinSeq + '\n')
    pro_f.close()
    cds.to_csv(path + specie + '.cds.csv', index=False, encoding='utf-8')
    
    
