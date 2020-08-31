import os, sys
import time
import warnings
import numpy as np
import pandas as pd
from sklearn.metrics import multilabel_confusion_matrix
from collections import Counter
import math
import itertools
from typing import *
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba, to_rgb
import seaborn as sns

rootdir = os.path.dirname(os.path.dirname(os.getcwd()))
sys.path.insert(0, rootdir)

from VCFPooling.poolSNPs import dataframe as vcfdf


# def ok_ko_mapper(i: str, labels: np.ndarray) -> dict:
#     d = {  # truth = i value
#         'okok': [],  # decoded True and imputed True
#         'okko': [],  # decoded True and imputed False  # should be 0!
#         'koko': [],  # decoded False and imputed False
#         'kook': [],  # decoded False and imputed True
#     }
#
#     d['okok'].append(str(i) + str(i) + str(i))
#
#     dko = list(labels)
#     dko.remove(str(i))
#     for j in dko:
#         d['okko'].append(str(i) + str(i) + str(j))
#         d['kook'].append(str(i) + str(j) + str(i))
#         for k in dko:
#             d['koko'].append(str(i) + str(j) + str(k))
#
#     return d
#
#
# def ok_ko_counter(dict_ok_ko: dict, counts: dict) -> dict:
#     '''counts = Counter(arr_str)'''
#     dcount = {}
#     for k, cases in dict_ok_ko.items():
#         dcount[k] = []
#         for cas in cases:
#             dcount[k].append(counts[cas])
#             dcount[cas] = counts[cas]
#     return dcount


# Params

imp_method = 'beagle'
genos = np.array([0, 1, 2, -1])
barcolors = ['#80013f',  # missing
             '#047495', '#00035b', '#748b97',  # full GT 0, 1, 2
             #'#dbb40c', '#c65102'  # half-missing GT
             ]
barcolors_rgb = [list(to_rgb(bc)) for bc in barcolors]
barcolors_rgba = []
for rgb in barcolors_rgb:
    barcolors_rgba.append(tuple(rgb + [0.5]))
    barcolors_rgba.append(tuple(rgb + [1.0]))

bin_step = 0.04
# x_bins = np.arange(0.0, 1.0 + bin_step, bin_step).round(decimals=2)
x_bins = [0.0, 0.02, 0.04, 0.06, 0.1, 0.2, 0.4, 0.6, 0.8, 0.9, 0.94, 0.96, 0.98, 1.0]
# lab_bins = np.arange(bin_step/2, 1.0, bin_step).round(decimals=2)
lab_bins = [0.01, 0.03, 0.05, 0.08, 0.15, 0.3, 0.5, 0.7, 0.85, 0.92, 0.95, 0.97, 0.99]

# Data

cube = np.random.choice(genos, size=(10000, 3), p=probs)

TP = np.array([[55, 0, 0, 5], [0, 0, 0, 30], [0, 0, 8, 2], [0, 0, 0, 0]])
PI = np.array([[55, 0, 0, 0], [0, 0, 0, 0], [0, 0, 8, 0], [20, 15, 2, 0]])

# true = '/home/camille/1000Genomes/src/VCFPooling/examples/IMP.chr20.snps.gt.vcf.gz'
# pooled = '/home/camille/1000Genomes/src/VCFPooling/examples/IMP.chr20.pooled.snps.gt.vcf.gz'
# imputed_beagle = '/home/camille/1000Genomes/src/VCFPooling/examples/IMP.chr20.pooled.imputed.vcf.gz'
true = '/home/camille/PoolImpHuman/data/20200722/IMP.chr20.snps.gt.vcf.gz'
pooled = '/home/camille/PoolImpHuman/data/20200812/IMP.chr20.pooled.snps.gt.vcf.gz'
imputed_beagle = '/home/camille/PoolImpHuman/data/20200722/IMP.chr20.pooled.imputed.vcf.gz'

dftrue = vcfdf.PandasMixedVCF(true, format='GT', indextype='chrom:pos')
dfpooled = vcfdf.PandasMixedVCF(pooled, format='GT', indextype='chrom:pos')
dfimputed= vcfdf.PandasMixedVCF(imputed_beagle, format='GT', indextype='chrom:pos')

af_bins = pd.cut(dftrue.af_info.values.squeeze(), bins=x_bins, labels=lab_bins)
binned_af = pd.Series(af_bins, index=dftrue.variants, name='binned_af')

dT = dftrue.trinary_encoding().join(binned_af)
dP = dfpooled.trinary_encoding().join(binned_af)
dI = dfimputed.trinary_encoding().join(binned_af)

dTgp = dT.groupby('binned_af')
dPgp = dP.groupby('binned_af')
dIgp = dI.groupby('binned_af')

genos2str = lambda x: ''.join(np.char.replace(x.astype(str), '-1', '-'))
genos_str = np.char.replace(genos.astype(str), '-1', '-')

arrdf = np.empty(
    (len(lab_bins) * (len(genos) ** 3), 5),
    dtype=int
)
dfcounts = pd.DataFrame(data=arrdf,
                        columns=['bin', 'trueGT', 'pooledGT', 'imputedGT', 'counts'],
                        dtype=int
                        )
dfcounts.bin = dfcounts.bin.astype(str)

b = 0
row = 0
for gT, gP, gI in zip(dTgp, dPgp, dIgp):  # bin by bin in the group
    xT = gT[1].drop('binned_af', axis=1).values.flatten()
    xP = gP[1].drop('binned_af', axis=1).values.flatten()
    xI = gI[1].drop('binned_af', axis=1).values.flatten()
    xA = np.array(list(zip(xT, xP, xI)))

    xA_str = np.apply_along_axis(genos2str, arr=xA, axis=-1)

    try:
        cntA = Counter(xA_str)
    except TypeError:  # if no markers in the AF-bin
        cntA = {}
    print(cntA)

    genos_counts = {}
    for i in genos_str:
        for j in genos_str:
            for k in genos_str:
                genos_counts[i + j + k] = 0
    genos_counts.update(cntA)

    for i in genos_str:
        for j in genos_str:
            for k in genos_str:
                dfcounts.loc[row, 'bin'] = str(lab_bins[b])
                dfcounts.loc[row, 'trueGT'] = -1 if i == '-' else int(i)
                dfcounts.loc[row, 'pooledGT'] = -1 if j == '-' else int(j)
                dfcounts.loc[row, 'imputedGT'] = -1 if k == '-' else int(k)
                dfcounts.loc[row, 'counts'] = genos_counts[i + j + k]
                row += 1
    b += 1

dfcounts = dfcounts.set_index('bin')
print(dfcounts)
dfcounts.to_csv('/home/camille/PoolImpHuman/confusion_true_pooled_imputed_beagle.csv', sep=",")
#TODO: save dfcounts to csv!


def get_dfbarplot_GT(gt: int):
    """data frame to bar plot"""
    dftrueGT = dfcounts[dfcounts.trueGT == gt][['pooledGT', 'imputedGT', 'counts']]
    decGT = dftrueGT[dftrueGT.pooledGT == gt][['imputedGT', 'counts']].groupby(['bin', 'imputedGT']).sum()
    nondecGT = dftrueGT[dftrueGT.pooledGT != gt][['imputedGT', 'counts']].groupby(['bin', 'imputedGT']).sum()
    decGTT = decGT.transpose()
    nondecGTT = nondecGT.transpose()
    dfplotGT = pd.concat([decGTT, nondecGTT]).reset_index(drop=True)
    dfplotGT['set'] = ['decoded', 'not decoded']
    dfplotGT.set_index('set', inplace=True)
    dfplotGT = dfplotGT.melt(value_name='counts',  ignore_index=False)
    dfplotGT.reset_index(inplace=True)
    dfplotGT = dfplotGT.pivot(index=['bin', 'set'], columns=['imputedGT'], values='counts')
    # dfplotGT_scaled = dfplotGT.divide(dfplotGT.sum(axis=1), axis=0)  # not, decoded sibe by side
    dfplotGT = dfplotGT.unstack()
    dfplotGT_scaled = dfplotGT.divide(dfplotGT.sum(axis=1), axis=0)

    return dfplotGT_scaled


dfplot0 = get_dfbarplot_GT(0)

plt.rcParams["figure.figsize"] = [12, 5]
fig, ax = plt.subplots()
dfplot0.plot(kind='bar', rot=45, stacked=True, color=barcolors_rgba, edgecolor='white')
plt.show()

plt.rcParams["figure.figsize"] = [4, 4]
for gx, g in enumerate(genos[:-1]):
    break
    dftrue = dfcounts[dfcounts.trueGT == g]
    sumcounts = sum(dftrue.counts)
    mincounts = min(dftrue.counts)
    maxcounts = max(dftrue.counts)
    #dftrue.counts = dftrue.counts / sumcounts
    ax = sns.scatterplot(data=dftrue, x='pooledGT', y='imputedGT',
                         #c=barcolors[gx], #  hue='counts',
                         size='counts', sizes=(mincounts, maxcounts))
    ax.set_xticks(genos)
    ax.set_yticks(genos)
    plt.title('trueGT = {}'.format(g))
    #plt.show()
