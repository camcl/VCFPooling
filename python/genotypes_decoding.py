"""
How does pooling decode the data? How "bad" does pooling make the data? What kind of missing data given the allele frequency?

Usage example:
$ python3 -u genotypes_decoding.py /home/camille/PoolImpHuman/data/20200812/IMP.chr20.snps.gt.vcf.gz /home/camille/PoolImpHuman/data/20200812/IMP.chr20.pooled.snps.gt.vcf.gz
"""

import numpy as np
import pandas as pd
from sklearn.metrics import multilabel_confusion_matrix, confusion_matrix
from sklearn.preprocessing import *
from scipy.stats import *
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, to_rgba_array, to_rgba
import seaborn as sns
import os, sys
import argparse

rootdir = os.path.dirname(os.path.dirname(os.getcwd()))
sys.path.insert(0, rootdir)

from VCFPooling.poolSNPs import dataframe as vcfdf


# Parse command-line arguments
parser = argparse.ArgumentParser(description='Barplots for genotypes states in pooled data')
parser.add_argument('truefile', metavar='truef', type=str, help='File with true genotypes GT', default=None)
parser.add_argument('pooledfile', metavar='pooledf', type=str, help='File with pooled genotypes decoded into GT', default=None)
argsin = parser.parse_args()

# Plotting features
true_genos = [0.0, 1.0, 2.0]
true_labels = ['0/0', '0/1', '1/1']
pooled_genos = [0.0, 1.0, 2.0, -0.5, 0.5, -1.0]
pooled_labels = ['0/0', '0/1', '1/1', '0/.', './1', './.']

# dashes_styles = [(0, ()), (0, ()), (0, ()),  # full GT
#                  (0, (5, 5)), (0, (5, 5)), (0, (1, 1))  # missing GT
#                  ]
dashes_styles = ['-', '-', '-',  # full GT
                 '--', '--', '-.'  # missing GT
                 ]
barcolors = ['#047495', '#00035b', '#748b97',
             '#dbb40c', '#c65102', '#80013f']
barcmap=ListedColormap([to_rgba(co) for co in barcolors])

bin_step = 0.04
# x_bins = np.arange(0.0, 1.0 + bin_step, bin_step).round(decimals=2)
x_bins = [0.0, 0.02, 0.04, 0.06, 0.1, 0.2, 0.4, 0.6, 0.8, 0.9, 0.94, 0.96, 0.98, 1.0]
# lab_bins = np.arange(bin_step/2, 1.0, bin_step).round(decimals=2)
lab_bins = [0.01, 0.03, 0.05, 0.08, 0.15, 0.3, 0.5, 0.7, 0.85, 0.92, 0.95, 0.97, 0.99]

truef = argsin.truefile
# argsin.truefile
# '/home/camille/1000Genomes/data/gt/stratified/IMP.chr20.snps.gt.chunk1000.vcf.gz'
# '/home/camille/1000Genomes/src/VCFPooling/examples/IMP.chr20.snps.gt.vcf.gz'
pooledf = argsin.pooledfile
# argsin.pooledfile
# '/home/camille/1000Genomes/data/gt/stratified/IMP.chr20.pooled.snps.gt.chunk1000.vcf.gz'
# '/home/camille/1000Genomes/src/VCFPooling/examples/IMP.chr20.pooled.snps.gt.vcf.gz'

dftrue = vcfdf.PandasMixedVCF(truef, format='GT')
dfpooled = vcfdf.PandasMixedVCF(pooledf, format='GT')
n_markers, n_samples = dftrue.genotypes().shape

# true = dftrue.hexa_encoding()
pooled = dfpooled.hexa_encoding()

af_bins = pd.cut(dftrue.aaf.values.squeeze(), bins=x_bins, labels=lab_bins)
binned_af = pd.Series(af_bins, index=dftrue.variants, name='binned_af')

# true = true.join(binned_af)
# true['dataset'] = ['true'] * true.shape[0]
# true = true.reset_index().set_index(['variants', 'binned_af', 'dataset'])

pooled = pooled.join(binned_af)
pooled['dataset'] = ['pooled'] * pooled.shape[0]
pooled = pooled.reset_index().set_index(['variants', 'binned_af', 'dataset'])

# Initialize counts for each AF bin and each genotype
dfcounts = pd.DataFrame(data=None, index=lab_bins, columns=pooled_genos)
for i in dfcounts.index:
    for j in dfcounts.columns:
        dfbins = pooled.loc[pooled.index.get_level_values('binned_af') == i]
        dfbins.reset_index(inplace=True, drop=True)
        counts_geno = dfbins.where(dfbins == j, axis=0).count()
        dfcounts.loc[i, j] = counts_geno.sum()
dfcounts.columns = pooled_labels


# scale result per bin
binscales = dfcounts.sum(axis=1)
dfcounts_scaled = pd.DataFrame(data=dfcounts.div(binscales, axis=0),  # minmax_scale(dfcounts.values, axis=0),
                               index=dfcounts.index,
                               columns=dfcounts.columns)

dfcounts_sized = dfcounts / (n_samples * n_markers)

figsize=4
plt.rcParams["figure.figsize"] = [figsize*3, figsize + 1]

ax = dfcounts_scaled.plot(kind='bar', stacked=True, rot=45,
                          color=barcolors, style=dashes_styles)  # cmap = sns.set_palette('GnBu_d')
ax.set_xlabel('Estimated alternate allele frequency')
ax.set_ylabel('Proportions of genotypes normalized per AF-bin')
plt.title('Genotypes proportions from pooled data in the study population')
plt.tight_layout()
plt.savefig(os.path.join(os.path.dirname(pooledf), 'genotypes_hexa_counts.pdf'))
plt.show()

ax_scaled = dfcounts_sized.plot(kind='bar', stacked=True, rot=45,
                                color=barcolors, style=dashes_styles)
ax_scaled.set_xlabel('Estimated alternate allele frequency')
ax_scaled.set_ylabel('Proportion of genotypes')
plt.title('Genotypes proportions from pooled data in the study population (total number of genotypes = {})'.format(n_samples * n_markers))
plt.tight_layout()
plt.savefig(os.path.join(os.path.dirname(pooledf), 'genotypes_hexa_proportions.pdf'))
plt.show()
