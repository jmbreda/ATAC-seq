import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import pyBigWig
import os

bin_size = 50
CHR=[str(c+1) for c in range(19)]+['X','Y']
infolder = 'results/mapping'

sample_metadata = pd.read_csv('resources/SraRunTable.txt')
# add a "Replicate" column
sample_metadata.loc[:,'replicate'] = [int(n.split('_')[-1][-1]) for n in sample_metadata.loc[:,'Library Name']]
Tissue = sample_metadata.TISSUE.unique()
Sex = ['Female','Male']

for tissue in Tissue:
    print(tissue)
    
    # load track tables
    fname = f'{infolder}/{tissue}_coverage.tsv'
    if os.path.exists(fname):
        f = open(fname, 'r')
        line1 = f.readline().replace('#','').replace("'",'').split()
        dtypes = {'chr': str, 'start': int, 'end': int}
        df = pd.read_csv(f, sep='\t', names=line1, dtype=dtypes,quotechar="'")

        # keep only chr in CHR
        idx_in = [c in CHR for c in df.chr]
        df = df.loc[idx_in,:]

        # sort table
        
        # make bigwigs header
        my_header = [(chr,df.loc[df.chr==chr,'end'].max()) for chr in CHR]
        
        # Make avg bigwigs
        for sex in Sex:
            print(f'{sex} average')
            if np.any( [sex in c for c in df.columns[3:]] ):
                f_out = f'{infolder}/{tissue}_{sex}_coverage.bw'
                with pyBigWig.open(f_out,'w') as bw:
                    bw.addHeader(my_header)
                    samp = [f'{sex}_{tissue}_Rep1',f'{sex}_{tissue}_Rep2']
                    for chr in CHR:
                        start=df.loc[df.chr==chr,'start'].values
                        end=df.loc[df.chr==chr,'end'].values
                        vals = df.loc[df.chr==chr,samp].values.mean(axis=1)
                        #sort
                        idx_sort = np.argsort(start)
                        bw.addEntries([chr]*len(start), start[idx_sort], ends=end[idx_sort], values=vals[idx_sort])

        # Make diff bigwig
        print('diff')
        if np.any( ['Female' in c for c in df.columns[3:]] ) & np.any( ['Male' in c for c in df.columns[3:]] ) :
            f_out = f'{infolder}/{tissue}_f_m_diff_coverage.bw'
            with pyBigWig.open(f_out,'w') as bw:
                bw.addHeader(my_header)
                samp_f = [f'Female_{tissue}_Rep1',f'Female_{tissue}_Rep2']
                samp_m = [f'Male_{tissue}_Rep1',f'Male_{tissue}_Rep2']
                for chr in CHR:
                    start=df.loc[df.chr==chr,'start'].values
                    end=df.loc[df.chr==chr,'end'].values
                    vals = df.loc[df.chr==chr,samp_f].values.mean(axis=1) - df.loc[df.chr==chr,samp_m].values.mean(axis=1)
                    #sort
                    idx_sort = np.argsort(start)
                    bw.addEntries([chr]*len(start), start[idx_sort], ends=end[idx_sort], values=vals[idx_sort])


