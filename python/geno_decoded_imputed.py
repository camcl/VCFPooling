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
parser = argparse.ArgumentParser(description='Barplots for genotypes states in imputed data')
parser.add_argument('truefile', metavar='truef', type=str, help='File with true genotypes GT', default=None)
parser.add_argument('pooledfile', metavar='pooledf', type=str, help='File with pooled genotypes decoded into GT', default=None)
parser.add_argument('imputedfile', metavar='imputedf', type=str, help='File with imputed genotypes decoded into GT', default=None)
parser.add_argument('outdir', metavar='outdir', type=str, help='Directory to save the plots', default=None)
parser.add_argument('algo', metavar='algo', type=str, help='Imputation mehtod', default=None)

argsin = parser.parse_args()

truef = argsin.truefile
pooledf = argsin.pooledfile
imputedf = argsin.imputedfile

outdir = argsin.outdir
if not os.path.exists(outdir):
    os.mkdir(outdir)
print('\r\nFigures will be saved in {}'.format(outdir).ljust(80, '.'))

algo = argsin.algo


# Plotting features and parameters
# Use trinary encoding even for pooled genotypes because the number of half-decoded genotypes is very low
true_genos = [0, 1, 2, -1]  # trinary encoding with int, not float
true_labels = ['0/0', '0/1', '1/1', './.']
pooled_genos = [0, 1, 2, -1]
pooled_labels = ['0/0', '0/1', '1/1', './.']
imputed_genos = [0, 1, 2, -1]
imputed_labels = ['0/0', '0/1', '1/1', './.']

# dashes_styles = [(0, ()), (0, ()), (0, ()),  # full GT
#                  (0, (5, 5)), (0, (5, 5)), (0, (1, 1))  # missing GT
#                  ]
dashes_styles = ['-', '-', '-',  # full GT
                 '--', '--', '-.'  # missing GT
                 ]
barcolors = ['#047495', '#00035b', '#748b97',  # full GT
             '#dbb40c', '#c65102', '#80013f'  # missing GT
             ]
barcmap = ListedColormap([to_rgba(co) for co in barcolors])

bin_step = 0.04
# x_bins = np.arange(0.0, 1.0 + bin_step, bin_step).round(decimals=2)
x_bins = [0.0, 0.02, 0.04, 0.06, 0.1, 0.2, 0.4, 0.6, 0.8, 0.9, 0.94, 0.96, 0.98, 1.0]
# lab_bins = np.arange(bin_step/2, 1.0, bin_step).round(decimals=2)
lab_bins = [0.01, 0.03, 0.05, 0.08, 0.15, 0.3, 0.5, 0.7, 0.85, 0.92, 0.95, 0.97, 0.99]

figsize=4
plt.rcParams["figure.figsize"] = [figsize*3, figsize + 1]
print('\r\nCounting genotypes'.ljust(80, '.'))


# Read and process data
print('\r\nReading data from {} and {}'.format(truef, imputedf).ljust(80, '.'))
dftrue = vcfdf.PandasMixedVCF(truef, format='GT', indextype='chrom:pos')
dfpooled = vcfdf.PandasMixedVCF(pooledf, format='GT', indextype='chrom:pos')
dfimputed = vcfdf.PandasMixedVCF(imputedf, format='GT', indextype='chrom:pos')
n_markers, n_samples = dftrue.genotypes().shape

# Use trinary encoding even for pooled genotypes because the number of half-decoded genotypes is very low
true = dftrue.trinary_encoding()
pooled = dfpooled.trinary_encoding()
imputed = dfimputed.trinary_encoding()

af_bins = pd.cut(dftrue.af_info.values.squeeze(), bins=x_bins, labels=lab_bins)
binned_af = pd.Series(af_bins, index=dftrue.variants, name='binned_af')

true = true.join(binned_af)
true['dataset'] = ['true'] * true.shape[0]
true = true.reset_index().set_index(['variants', 'binned_af', 'dataset'])

pooled = pooled.join(binned_af)
pooled['dataset'] = ['pooled'] * pooled.shape[0]
pooled = pooled.reset_index().set_index(['variants', 'binned_af', 'dataset'])

imputed = imputed.join(binned_af)
imputed['dataset'] = ['imputed'] * imputed.shape[0]
imputed = imputed.reset_index().set_index(['variants', 'binned_af', 'dataset'])

# Initialize counts for each AF bin and each genotype: form confusion matrices
print('\r\nCounting genotypes'.ljust(80, '.'))
cmGT00 = np.zeros((len(lab_bins), len(imputed_genos)), dtype=float)  # (13, 3) shape
cmGT01 = np.zeros((len(lab_bins), len(imputed_genos)), dtype=float)  # (13, 3) shape
cmGT11 = np.zeros((len(lab_bins), len(imputed_genos)), dtype=float)  # (13, 3) shape

cmatrix = np.zeros((
    len(lab_bins),  # dim0
    len(true_genos),  # dim1
    len(pooled_genos),  # dim2
    len(imputed_genos)  # dim3
), dtype=float)

for bx, bin in enumerate(lab_bins):
        #bin = lab_bins[i]
        dT = true.loc[true.index.get_level_values('binned_af') == bin]
        dP = pooled.loc[pooled.index.get_level_values('binned_af') == bin]
        dI = imputed.loc[imputed.index.get_level_values('binned_af') == bin]

        print(dT)
        print('\n')
        print(dP)
        print('\n')
        print(dI)

        confusionTP = np.zeros((len(true_genos), len(pooled_genos)))
        for y_true, y_pred in zip(dT.values, dP.values):
            cm = confusion_matrix(y_true.astype(str), y_pred.astype(str),
                                  labels=np.array(pooled_genos).astype(str))
            confusionTP = confusionTP + cm

        confusionPI = np.zeros((len(pooled_genos), len(imputed_genos)))
        for y_true, y_pred in zip(dP.values, dI.values):
            cm = confusion_matrix(y_true.astype(str), y_pred.astype(str),
                                  labels=np.array(imputed_genos).astype(str))
            confusionPI = confusionPI + cm

        for ix, i in enumerate(true_genos):
            for jx, j in enumerate(pooled_genos):
                for kx, k in enumerate(imputed_genos):
                    cmatrix[bx, ix, jx, kx] = 0.0

        #
        # cross_genos = confusion[-3:]  # rectangle, no missing genotypes in the true dataset
        # scaled_genos = minmax_scale(cross_genos, axis=1)
        #
        # cmGT00[i] = cross_genos[0, :] / sum(cross_genos[0, :])
        # cmGT01[i] = cross_genos[1, :] / sum(cross_genos[1, :])
        # cmGT11[i] = cross_genos[2, :] / sum(cross_genos[2, :])

dfGT00 = pd.DataFrame(data=cmGT00, index=lab_bins, columns=imputed_labels)
# dfGT00['dataset'] = 'beagle'
dfGT01 = pd.DataFrame(data=cmGT01, index=lab_bins, columns=imputed_labels)
# dfGT01['dataset'] = 'beagle'
dfGT11 = pd.DataFrame(data=cmGT11, index=lab_bins, columns=imputed_labels)
# dfGT11['dataset'] = 'beagle'

print('\n')
print(dfGT00)
print('\n')
print(dfGT01)
print('\n')
print(dfGT11)
print('\n')


# Plot processed data
print('\r\nPlotting results'.ljust(80, '.'))
ax00 = dfGT00.plot(kind='bar', stacked=True, rot=45,
                   color=barcolors, style=dashes_styles)
ax00.set_xlabel('True alternate allele frequency')
ax00.set_ylabel('Proportions of genotypes scaled per AAF-bin')
plt.title('Genotypes imputed for true 0/0 genotype in the study population')
plt.tight_layout()
plt.savefig(os.path.join(outdir, 'true00_geno_decode_imputed_scaled_proportions.pdf'))
plt.show()

ax01 = dfGT01.plot(kind='bar', stacked=True, rot=45,
                   color=barcolors, style=dashes_styles)
ax01.set_xlabel('True alternate allele frequency')
ax01.set_ylabel('Proportions of genotypes scaled per AAF-bin')
plt.title('Genotypes imputed for true 0/1 genotype in the study population')
plt.tight_layout()
plt.savefig(os.path.join(outdir, 'true01_geno_decode_imputed_scaled_proportions.pdf'))
plt.show()

ax11 = dfGT11.plot(kind='bar', stacked=True, rot=45,
                   color=barcolors, style=dashes_styles)
ax11.set_xlabel('True alternate allele frequency')
ax11.set_ylabel('Proportions of genotypes scaled per AAF-bin')
plt.title('Genotypes imputed for true 1/1 genotype in the study population')
plt.tight_layout()
plt.savefig(os.path.join(outdir, 'true11_geno_decode_imputed_scaled_proportions.pdf'))
plt.show()


# All plots in same page
n_rows = len(imputed_genos)
n_cols = 1
fig, axes = plt.subplots(n_rows, n_cols, figsize=[figsize*3*n_cols, (figsize + 1)*n_rows], constrained_layout=True)

_ = dfGT00.plot(kind='bar', stacked=True, rot=45, ax=axes[0],
                color=barcolors, style=dashes_styles)
axes[0].set_xlabel('True alternate allele frequency')
axes[0].set_ylabel('Proportions of genotypes scaled per AAF-bin')
axes[0].set_title('Genotypes imputed for true 0/0 genotype in the study population')

_ = dfGT01.plot(kind='bar', stacked=True, rot=45, ax=axes[1],
                color=barcolors, style=dashes_styles)
axes[1].set_xlabel('True alternate allele frequency')
axes[1].set_ylabel('Proportions of genotypes scaled per AAF-bin')
axes[1].set_title('Genotypes imputed for true 0/1 genotype in the study population')

_ = dfGT11.plot(kind='bar', stacked=True, rot=45, ax=axes[2],
                color=barcolors, style=dashes_styles)
axes[2].set_xlabel('True alternate allele frequency')
axes[2].set_ylabel('Proportions of genotypes scaled per AAF-bin')
axes[2].set_title('Genotypes imputed for true 1/1 genotype in the study population')

fig.suptitle(algo)
# plt.subplots_adjust(top=0.85)
fig.savefig(os.path.join(outdir, '{}_geno_decode_imputed_scaled_proportions.pdf'.format(algo)))
