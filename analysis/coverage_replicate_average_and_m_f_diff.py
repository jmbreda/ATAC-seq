import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import pyBigWig

bin_size = 50
CHR=[str(c+1) for c in range(19)]+['X','Y']
infolder = 'results/mapping'
Samples = {'liver_male_1':'SRR8119839_coverage.bw',
           'liver_male_2':'SRR8119838_coverage.bw',
           'liver_female_1':'SRR8119852_coverage.bw',
           'liver_female_2':'SRR8119853_coverage.bw'}

# TODO 
# intersect bw_rep1.interval(chr) and bw_rep2.interval(chr) and get average coverage.
# get scatter and diff bw_male_avg.interval(chr) bw_female_avg.interval(chr)

samples = pd.read_csv('resources/SraRunTable.txt')
# add a "Replicate" column
samples.loc[:,'replicate'] = [int(n.split('_')[-1][-1]) for n in samples.loc[:,'Library Name']]
# filter for my_samples
my_samples = ['SRR8119838','SRR8119839','SRR8119852','SRR8119853'] 
idx = [s in my_samples for s in samples.Run]
samples = samples.iloc[idx,:]
Tissue = samples.TISSUE.unique()
Sex = samples.sex.unique()

for tissue in Tissue:
    print(tissue)
    # initialize interval starts-ends
    start = {}
    end = {}
    for chr in CHR:
        start[chr] = set()
        end[chr] = set()

    bw = {}
    for sex in Sex:
        print(sex)
        for r in [1,2]:
            print(r)
            f = samples.loc[(samples.TISSUE==tissue) & (samples.sex == sex) & (samples.replicate==r),'Run'].values
            if len(f) < 1:
                print(f'sample ({tissue} {sex} {r}) not found')
            elif len(f) > 2:
                print(f'multiple samples ({tissue} {sex} {r}) found')
            else:
                # load bw file
                bw[f'{sex}_{r}'] = pyBigWig.open(f'{infolder}/{f[0]}_coverage.bw')

                #now add intarvals
                for chr in CHR:
                    print(chr)
                    start[chr].update( [i[0] for i in bw[f'{sex}_{r}'].intervals(chr)] )
                    end[chr].update( [i[1] for i in bw[f'{sex}_{r}'].intervals(chr)] )
    
    # sort starts and ends
    for chr in CHR:
        start[chr] = sorted(start[chr])
        end[chr] = sorted(end[chr])

    # make header with male_1 length
    my_header = [(chr,bw['male_1'].chroms()[chr]) for chr in CHR]
    
    # Make avg bigwigs
    for sex in Sex:
        f_out = f'{infolder}/{tissue}_{sex}_coverage.bw'
        bw[sex] = pyBigWig.open(f_out,'w')
        bw[sex].addHeader(my_header)
        bw[sex].addheader()
        for chr in CHR:
            vals = .5*( bw[f'{sex}_1'].intervals(chr) + bw[f'{sex}_2'].intervals(chr) )
            bw.addEntries(chr, start[chr], ends=end[chr], values=vals)

    # Make diff bigwig

