#sjd_sle_overlaps.py

import pandas as pd
from intervaltree import IntervalTree
from intervaltree import Interval
import numpy as np

chr = 'CHR'
pos = 'POS'
offset = 250000

sjd_table = pd.read_csv('../Tables/SjD_23.csv', lineterminator='\n')
sle_table = pd.read_csv('../Tables/SLE_182.csv', lineterminator='\n')

#rename columns CHR and POS and drop the \r in column name
sjd_table.columns = ['Variant and risk allele', 'P-value', 'P-value annotation', 'RAF', 'OR', 'Beta', 'CI', 'Mapped gene', 'Reported trait', 'Trait(s)','Background trait(s)', 'Study accession', 'PMID', 'Unnamed: 13','Location', 'CHR', 'POS', 'Putative Causal Gene','supporting evidence', 'source', 'P']
#replace the carriage returns with null string
sjd_table['P'] = sjd_table['P'].str.replace(r'\r', '')

#rename columns CHR and POS and drop the \r in column name
sle_table.columns = ['Variant and risk allele', 'P-value', 'P-value annotation', 'RAF', 'OR','Beta', 'CI', 'Mapped gene', 'Reported trait', 'Trait(s)','Background trait(s)', 'Study accession', 'PubMed ID', 'First Author','Location', 'P', 'CHR', 'POS', 'Region (Value) ','Final Table?', 'Putative Causal Gene', 'OpenTargets']
#replace the carriage returns with null string
sle_table['OpenTargets'] = sle_table['OpenTargets'].str.replace(r'\r', '')
#rearrange columns the URL at the end is unsightly.
sle_table = sle_table[['Variant and risk allele', 'P-value', 'P-value annotation', 'RAF', 'OR','Beta', 'CI', 'Mapped gene', 'Reported trait', 'Trait(s)','Background trait(s)', 'Study accession', 'PubMed ID', 'First Author','Location', 'P', 'CHR', 'POS',  'OpenTargets', 'Region (Value) ','Final Table?', 'Putative Causal Gene']]

remap_dict = {'X':23}
sjd_table[chr] = sjd_table[chr].replace(remap_dict)
sjd_table[chr] = sjd_table[chr].astype('int64')
sjd_table[pos] = sjd_table[pos].astype('int64')

sle_table[chr] = sle_table[chr].replace(remap_dict)
sle_table[chr] = sle_table[chr].astype('int64')
sle_table[pos] = sle_table[pos].astype('int64')

sjd_table = sjd_table.sort_values(['CHR','POS'])
sle_table = sle_table.sort_values(['CHR','POS'])


#build sle interval tree [build the bigger tree – building should be cheaper than querying]
sledf = sle_table.groupby(chr)
sle_chromosome_list = sle_table.CHR.unique()
sle_trees = [IntervalTree() for i in range(len(sle_chromosome_list)+1)]

for key, sle_chr in sledf:
    print(int(key))
    for row in sle_chr.itertuples():
        region_start = row.POS-offset
        region_end = row.POS+offset
        #each entry is POS start, POS end, chr
        #key is offset by 1 chr 1-23, array 0 to 22
        sle_trees[int(key)-1].addi(region_start,region_end)

sjddf = sjd_table.groupby(chr)
sjd_chromosome_list = sjd_table.CHR.unique()
sjd_trees = [IntervalTree() for i in range(len(sjd_chromosome_list)+1)]

overlap_count = 0
for key, sjd_chr in sjddf:
    for row in sjd_chr.itertuples():
        region_start = row.POS-offset
        region_end = row.POS+offset
        if sle_trees[int(key)-1].overlap(region_start,region_end):
            overlap_count = overlap_count + 1
            print("overlap: chr, interval ")
            print(int(key)-1)
            print(sle_trees[int(key)-1].overlap(region_start,region_end))

print("overlapping region #, k:")
print(overlap_count)
print("sjogren region #, s:")
r, c = sjd_table.shape
print(r)
r, c = sle_table.shape
print("sle region #, M:")
print(r)

N = np.floor((3298912062/500000)+1)
print("N: ")
print (N)
'''
k = 18
s = 23
M = 182
N = 6598

Parameters: 18, 23, 182, 6598
expected number of successes = 0.634434677174902
the results are over enriched 28.37 fold compared to expectations
hypergeometric p-value = 1.095536032654934e-24
'''
