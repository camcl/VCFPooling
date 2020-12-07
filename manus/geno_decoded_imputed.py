"""
What kind of genotype confusion/cooncordance/discordance per genotype state given the allelic frequency?
Between true and pooled and imputed.
(not minor allele because different kind of homozygotes)

Usage example:
(Beagle)
$ python3 -u geno_decode_imputed.py /home/camille/PoolImpHuman/data/20200722/IMP.chr20.snps.gt.vcf.gz /home/camille/PoolImpHuman/data/20200812/IMP.chr20.pooled.snps.gt.vcf.gz /home/camille/PoolImpHuman/data/20200722/IMP.chr20.pooled.imputed.vcf.gz /home/camille/PoolImpHuman/results/20200722 beagle

(Phaser)
$ python3 -u geno_decode_imputed.py /home/camille/PoolImpHuman/data/20200817/IMP.chr20.snps.gt.vcf.gz /home/camille/PoolImpHuman/data/20200812/IMP.chr20.pooled.snps.gt.vcf.gz /home/camille/PoolImpHuman/data/20200817/IMP.chr20.pooled.imputed.vcf.gz /home/camille/PoolImpHuman/results/20200817 phaser

(Beagle example subset)
VCFPooling/examples$ python3 -u ../python/geno_decode_imputed.py ./IMP.chr20.snps.gt.vcf.gz ./IMP.chr20.pooled.snps.gt.vcf.gz ./IMP.chr20.pooled.imputed.vcf.gz /.
"""

import numpy as np
import pandas as pd
from sklearn.metrics import multilabel_confusion_matrix, confusion_matrix
from sklearn.preprocessing import *
from collections import Counter
from scipy.stats import *
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, to_rgb, to_rgba
import seaborn as sns
import os, sys
import argparse

rootdir = os.path.dirname(os.path.dirname(os.getcwd()))
sys.path.insert(0, rootdir)

from VCFPooling.poolSNPs import dataframe as vcfdf


# Parse command-line arguments
# parser = argparse.ArgumentParser(description='Barplots for genotypes states in imputed data')
# parser.add_argument('truefile', metavar='truef', type=str, help='File with true genotypes GT', default=None)
# parser.add_argument('pooledfile', metavar='pooledf', type=str, help='File with pooled genotypes decoded into GT', default=None)
# parser.add_argument('imputedfile', metavar='imputedf', type=str, help='File with imputed genotypes decoded into GT', default=None)
# parser.add_argument('outdir', metavar='outdir', type=str, help='Directory to save the plots', default=None)
# parser.add_argument('algo', metavar='algo', type=str, help='Imputation mehtod', default=None)
#
# argsin = parser.parse_args()
#
# truef = argsin.truefile
# pooledf = argsin.pooledfile
# imputedf = argsin.imputedfile
#
# outdir = argsin.outdir
# if not os.path.exists(outdir):
#     os.mkdir(outdir)
# print('\r\nFigures will be saved in {}'.format(outdir).ljust(80, '.'))
#
# algo = argsin.algo


# Plotting features and parameters

# Use trinary encoding even for pooled genotypes because the number of half-decoded genotypes is very low
genos = np.array([0, 1, 2, -1])  # trinary encoding with int, not float
genos_labels = ['0/0', '0/1', '1/1', './.']
genos_int2label = dict(zip(genos, genos_labels))
genos_label2int = dict(zip(genos_labels, genos))
genos2str = lambda x: ''.join(np.char.replace(x.astype(str), '-1', '-'))
genos_str = np.char.replace(genos.astype(str), '-1', '-')

# algo = 'phaser'
algo = 'beagle'
table_counts = 'confusion_true_pooled_imputed_{}.csv'.format(algo)
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

# bin_step = 0.04
# x_bins = np.arange(0.0, 1.0 + bin_step, bin_step).round(decimals=2)
x_bins = [0.0, 0.02, 0.04, 0.06, 0.1, 0.2, 0.4, 0.6, 0.8, 0.9, 0.94, 0.96, 0.98, 1.0]
# lab_bins = np.arange(bin_step/2, 1.0, bin_step).round(decimals=2)
lab_bins = [0.01, 0.03, 0.05, 0.08, 0.15, 0.3, 0.5, 0.7, 0.85, 0.92, 0.95, 0.97, 0.99]

figsize=4
plt.rcParams["figure.figsize"] = [figsize*3, figsize + 1]


# Read and process data

# print('\r\nReading data from {} and {}'.format(truef, imputedf).ljust(80, '.'))
# dftrue = vcfdf.PandasMixedVCF(truef, format='GT', indextype='chrom:pos')
# dfpooled = vcfdf.PandasMixedVCF(pooledf, format='GT', indextype='chrom:pos')
# dfimputed = vcfdf.PandasMixedVCF(imputedf, format='GT', indextype='chrom:pos')
# n_markers, n_samples = dftrue.genotypes().shape

# Use trinary encoding even for pooled genotypes because the number of half-decoded genotypes is very low
# true = dftrue.trinary_encoding()
# pooled = dfpooled.trinary_encoding()
# imputed = dfimputed.trinary_encoding()

truef = '/home/camille/PoolImpHuman/data/20200722/IMP.chr20.snps.gt.vcf.gz'
pooledf = '/home/camille/PoolImpHuman/data/20200812/IMP.chr20.pooled.snps.gt.vcf.gz'
# imputedf = '/home/camille/PoolImpHuman/data/20200817/IMP.chr20.pooled.imputed.vcf.gz'  # phaser
imputedf = '/home/camille/PoolImpHuman/data/20200722/IMP.chr20.pooled.imputed.vcf.gz'  # beagle
# outdir = '/home/camille/PoolImpHuman/results/20200817'  # phaser
outdir = '/home/camille/PoolImpHuman/results/20200722'  # beagle

dftrue = vcfdf.PandasMixedVCF(truef, format='GT', indextype='chrom:pos')
dfpooled = vcfdf.PandasMixedVCF(pooledf, format='GT', indextype='chrom:pos')
dfimputed= vcfdf.PandasMixedVCF(imputedf, format='GT', indextype='chrom:pos')

af_bins = pd.cut(dftrue.af_info.values.squeeze(), bins=x_bins, labels=lab_bins)
binned_af = pd.Series(af_bins, index=dftrue.variants, name='binned_af')

if not os.path.exists(os.path.join(outdir, table_counts)):  # no counts values written in csv file yet
    dT = dftrue.trinary_encoding().join(binned_af)
    dP = dfpooled.trinary_encoding().join(binned_af)
    dI = dfimputed.trinary_encoding().join(binned_af)

    dTgp = dT.groupby('binned_af')
    dPgp = dP.groupby('binned_af')
    dIgp = dI.groupby('binned_af')


    # Initialize counts for each AF bin and each genotype

    print('\r\nCounting genotypes'.ljust(80, '.'))
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
    dfcounts.to_csv(os.path.join(outdir, table_counts), sep=",")

else:  # counts already in table, save time
    dfcounts = pd.read_csv(os.path.join(outdir, table_counts), sep=',')

dfcounts.reset_index(inplace=True)
print('Genotypes counts:\r\n', dfcounts)


# Plot processed data

print('\r\nPlotting results'.ljust(80, '.'))


def get_dfbarplot_GT(df: pd.DataFrame, gt: int) -> pd.DataFrame:
    """data frame for bar plot"""
    print(df)
    dftrueGT = df[df.trueGT == gt][['bin', 'pooledGT', 'imputedGT', 'counts']]
    decGT = dftrueGT[dftrueGT.pooledGT == gt][['bin', 'imputedGT', 'counts']].groupby(['bin', 'imputedGT']).sum()
    nondecGT = dftrueGT[dftrueGT.pooledGT != gt][['bin', 'imputedGT', 'counts']].groupby(['bin', 'imputedGT']).sum()

    dfplotGT = pd.concat([decGT.transpose(), nondecGT.transpose()]).reset_index(drop=True)
    dfplotGT['set'] = ['decoded', 'not decoded']
    dfplotGT.set_index('set', inplace=True)
    dfplotGT = dfplotGT.melt(value_name='counts',  ignore_index=False)
    dfplotGT.reset_index(inplace=True)
    dfplotGT = dfplotGT.pivot(index=['bin', 'set'], columns=['imputedGT'], values='counts')
    # dfplotGT_scaled = dfplotGT.divide(dfplotGT.sum(axis=1), axis=0)  # not, decoded sibe by side
    dfplotGT = dfplotGT.unstack()
    dfplotGT_scaled = dfplotGT.divide(dfplotGT.sum(axis=1), axis=0)

    return dfplotGT_scaled



# Single figure per true genotype
for gx, g in enumerate(genos[:-1]):
    dfplot = get_dfbarplot_GT(dfcounts, gx)

    fig, ax = plt.subplots()
    _ = dfplot.plot(kind='bar', ax=ax, rot=45, stacked=True, color=barcolors_rgba)  #, edgecolor='white')
    ax.set_title('Pooling and imputation outcomes for true {} genotype in the study population'.format(genos_labels[gx]))
    ax.set_xlabel('True alternate allele frequency')
    ax.set_ylabel('Proportions of genotypes scaled per AAF-bin')
    handles, labels = ax.get_legend_handles_labels()
    new_labels = []
    for leg_lab in labels:
        a = np.array(leg_lab.strip('()').split(','))
        new_labels.append('{} to {} and imputed to {}'.format(a[-1], genos_labels[gx], genos_int2label[int(a[0])]))
        # ex. '(2, not decoded)' -> '0 not decoded and imputed to 2' with gx=0
    ax.legend(bbox_to_anchor=(1.05, 0.5), loc='center left')
    ax.legend(handles, new_labels, bbox_to_anchor=(1.05, 0.5), loc='center left')
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, 'true{}_geno_decoded_imputed_confusion.pdf'.format(g)))
    plt.show()

# All plots in same page
n_rows = len(genos[:-1])
n_cols = 1
fig, axes = plt.subplots(n_rows, n_cols, figsize=[figsize*3*n_cols, (figsize + 1)*n_rows], constrained_layout=True)

for gx, g in enumerate(genos[:-1]):
    dfplot = get_dfbarplot_GT(dfcounts, gx)

    ax = axes[gx]
    _ = dfplot.plot(kind='bar', ax=ax, rot=45, stacked=True, color=barcolors_rgba)  #, edgecolor='white')
    ax.set_title('Pooling and imputation outcomes for true {} genotype in the study population'.format(genos_labels[gx]))
    ax.set_xlabel('True alternate allele frequency')
    ax.set_ylabel('Proportions of genotypes scaled per AAF-bin')
    box = ax.get_position()
    ax.set_position([box.x0 * 0.85, box.y0, box.width * 0.7, box.height * 0.8])  # resize position
    handles, labels = ax.get_legend_handles_labels()
    new_labels = []
    for leg_lab in labels:
        a = np.array(leg_lab.strip('()').split(','))
        new_labels.append('{} to {} and imputed to {}'.format(a[-1], genos_labels[gx], genos_int2label[int(a[0])]))
        # ex. '(2, not decoded)' -> '0 not decoded and imputed to 2' with gx=0
    ax.legend(bbox_to_anchor=(1.05, 0.5), loc='center left')
    ax.legend(handles, new_labels, bbox_to_anchor=(1.05, 0.5), loc='center left')
plt.suptitle(algo)  # NOT fig.suptitle(): it collapses the title with the plots
fig.savefig(os.path.join(outdir, '{}_geno_decoded_imputed_confusion.pdf'.format(algo)))
plt.show()


# # Extra tools
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