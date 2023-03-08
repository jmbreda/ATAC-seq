import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import pyBigWig

bin_size = 50
CHR=[str(c+1) for c in range(19)]+['X','Y']
infolder = 'results/mapping'

sample_metadata = pd.read_csv('resources/SraRunTable.txt')
# add a "Replicate" column
samples.loc[:,'replicate'] = [int(n.split('_')[-1][-1]) for n in samples.loc[:,'Library Name']]
# filter for my_samples
#my_samples = ['SRR8119838','SRR8119839','SRR8119852','SRR8119853']
#idx = [s in my_samples for s in samples.Run]
samples = samples.iloc[idx,:]
Tissue = samples.TISSUE.unique()
Sex = samples.sex.unique()
Tissue = ['Liver']

for tissue in Tissue:
    print(tissue)
    
    # load track tables
    dtypes = {chr: srt, start: int, end: int}
    fname = f'{infolder}/{tissue}_coverage.tsv'
    f = open(fname, 'r')
    line1 = f.readline().replace('#', '').replace("'",'').split()
    df = pd.read_csv(f, sep='\t', names=line1, dtype=dtypes,quotechar="'")
    
    # make header with male_1 length
    #my_header = [(chr,df[df.chr==chr,'end']) for chr in CHR]
    
    # Make avg bigwigs
    #for sex in Sex:
        #f_out = f'{infolder}/{tissue}_{sex}_coverage.bw'
        #bw[sex] = pyBigWig.open(f_out,'w')
        #bw[sex].addHeader(my_header)
        #bw[sex].addheader()
        #for chr in CHR:
        #    vals = bw[f'{sex}_1'].intervals(chr)
        #    bw.addEntries(chr, start[chr], ends=end[chr], values=)

    # Make diff bigwig

