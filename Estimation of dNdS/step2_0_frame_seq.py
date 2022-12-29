import pandas as pd
import re
path = '/media/hp/disk4/YT/protein_homology/ref/CRIF/'
species = ['hg38', 'mus']
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
        'UGC': 'C', 'UGU': 'C', 'UGA': '*', 'UGG': 'W',
        }
for specie in species:
    pep_in = pd.read_csv(path + specie + '_mix.csv', encoding='utf-8')
    p_out = open(path + specie + '_p0_protein.fa', 'w')
    m_out = open(path + specie + '_m0_protein.fa', 'w')
    ID = pep_in['id']
    start = pep_in['start']
    stop = pep_in['stop']
    mstop = pep_in[' +1_stop']
    pstop = pep_in['-1_stop']
    SEQ = pep_in['SEQ']
    motif = pep_in['motif']
    
    ID_v = list(map(str, ID))
    start_v = list(map(int, start))
    stop_v = list(map(int, stop))
    pstop_v = list(map(int, pstop))
    mstop_v = list(map(int, mstop))
    SEQ_v = list(map(str, SEQ))
    motif_v = list(map(str, motif))
    
    for i in range(0, len(ID_v)):
        p_rna = re.sub('T', 'U', SEQ_v[i][(stop_v[i] - 9) : pstop_v[i]])
        m_rna = re.sub('T', 'U', SEQ_v[i][(stop_v[i] - 9) : mstop_v[i]])
        proteinSeq = ""
        for codonStart in range(0, len(p_rna), 3):
            codon = p_rna[codonStart:codonStart + 3]
            if codon in codonTable:
                proteinSeq += codonTable[codon]
        p_out.write('>' + ID_v[i] + '+' + motif_v[i] + '+' + str(start_v[i]) + '+' + str(stop_v[i]) + '\n' + proteinSeq + '\n')

        for codonStart in range(0, len(m_rna), 3):
            codon = m_rna[codonStart:codonStart + 3]
            if codon in codonTable:
                proteinSeq += codonTable[codon]
        m_out.write('>' + ID_v[i] + '+' + motif_v[i] + '+' + str(start_v[i]) + '+' + str(stop_v[i]) + '\n' + proteinSeq + '\n')   
    m_out.close()
    p_out.close()    
