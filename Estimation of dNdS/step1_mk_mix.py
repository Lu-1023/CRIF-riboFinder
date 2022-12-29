import pandas as pd
path = '/media/hp/disk4/YT/protein_homology/ref/CRIF/'
species = ['hg38', 'mus']
for specie in species:
    pep_txt =pd.read_csv(path + 'pep_' + specie + '.txt', delimiter='\t', header=0, encoding='utf-8')
    pep_txt = pep_txt.drop_duplicates()
    cds_csv = pd.read_csv(path + specie + '.cds.csv', encoding='utf-8')
    cds_csv = cds_csv.drop_duplicates()
    mix =  pd.merge(pep_txt, cds_csv, how='left', on='id')
    mix.to_csv(path + specie + '_mix.csv', index=False, encoding='utf-8')
