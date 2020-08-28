"""
How does imputation confuse the data?

Usage example (from Beagle's imputed genotypes):
$ python3 -u imputation_confusion.py /home/camille/PoolImpHuman/data/20200722/IMP.chr20.snps.gt.vcf.gz /home/camille/PoolImpHuman/data/20200722/IMP.chr20.pooled.imputed.vcf.gz /home/camille/PoolImpHuman/results/20200722

Usage example (from Phaser's imputed genotypes):
$ python3 -u imputation_confusion.py /home/camille/PoolImpHuman/data/20200722/IMP.chr20.snps.gt.vcf.gz /home/camille/PoolImpHuman/data/20200817/IMP.chr20.pooled.imputed.vcf.gz /home/camille/PoolImpHuman/results/20200817
"""

import numpy as np
import pandas as pd
from sklearn.metrics import multilabel_confusion_matrix, confusion_matrix
from sklearn.preprocessing import *
from scipy.stats import *
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import seaborn as sns
import os, sys
import argparse

rootdir = os.path.dirname(os.path.dirname(os.getcwd()))
sys.path.insert(0, rootdir)

from VCFPooling.poolSNPs import dataframe as vcfdf


# Parse command-line arguments
parser = argparse.ArgumentParser(description='Confusion matrices for genotypes states in imputed data')
parser.add_argument('truefile', metavar='truef', type=str, help='File with true genotypes GT', default=None)
parser.add_argument('imputedfile', metavar='imputedf', type=str, help='File with imputed genotypes in GT format', default=None)
parser.add_argument('outdir', metavar='outdir', type=str, help='Directory to save the plots', default=None)
argsin = parser.parse_args()

truef = argsin.truefile
imputedf = argsin.imputedfile
outdir = argsin.outdir
if not os.path.exists(outdir):
    os.mkdir(outdir)
print('\r\nFigures will be saved in {}'.format(outdir).ljust(80, '.'))


# Plotting features and parameters
true_genos = [0.0, 1.0, 2.0]
true_labels = ['0/0', '0/1', '1/1']
imputed_genos = [-1.0, -0.5, 0.5, 0.0, 1.0, 2.0]
imputed_labels = ['./.', '0/.', './1', '0/0', '0/1', '1/1']

heatcolors = sns.cubehelix_palette(n_colors=10)
heatmask = np.array([[True, True, True, False, False, False],
                     [True, True, True, False, False, False],
                     [True, True, True, False, False, False]])  # mask impossible decoding combinations in cross_genos

x_bins = [0.0, 0.02, 0.04, 0.06, 0.1, 0.2, 0.4, 0.6, 0.8, 0.9, 0.94, 0.96, 0.98, 1.0]
lab_bins = [0.01, 0.03, 0.05, 0.08, 0.15, 0.3, 0.5, 0.7, 0.85, 0.92, 0.95, 0.97, 0.99]

figsize=4
plt.rcParams["figure.figsize"] = [figsize*(len(x_bins)//2)*2, figsize*2]


# Read and process data
print('\r\nReading data from {} and {}'.format(truef, imputedf).ljust(80, '.'))
dftrue = vcfdf.PandasMixedVCF(truef, format='GT', indextype='chrom:pos')
dfimputed = vcfdf.PandasMixedVCF(imputedf, format='GT', indextype='chrom:pos')
n_markers, n_samples = dftrue.genotypes().shape

true = dftrue.hexa_encoding()
imputed = dfimputed.hexa_encoding()

af_bins = pd.cut(dftrue.af_info.values.squeeze(), bins=x_bins, labels=lab_bins)
binned_af = pd.Series(af_bins, index=dftrue.variants, name='binned_af')

true = true.join(binned_af)
true['dataset'] = ['true'] * true.shape[0]
true = true.reset_index().set_index(['variants', 'binned_af', 'dataset'])

imputed = imputed.join(binned_af)
imputed['dataset'] = ['imputed'] * imputed.shape[0]
imputed = imputed.reset_index().set_index(['variants', 'binned_af', 'dataset'])


# Plot processed data
print('\r\nCounting genotypes and plotting results'.ljust(80, '.'))
fig, axes = plt.subplots(2, len(x_bins)//2)

for i, ax in enumerate(axes.flatten()):
    if i < len(lab_bins):
        bin = lab_bins[i]
        dT = true.loc[true.index.get_level_values('binned_af') == bin]
        dP = imputed.loc[imputed.index.get_level_values('binned_af') == bin]

        confusion = np.zeros((len(imputed_genos), len(imputed_genos)))
        for y_true, y_pred in zip(dT.values, dP.values):
            cm = confusion_matrix(y_true.astype(str), y_pred.astype(str),
                                  labels=np.array(imputed_genos).astype(str))
            confusion = confusion + cm

        cross_genos = confusion[-3:]  # rectanle, no missing genotypes in the true dataset
        scaled_genos = minmax_scale(cross_genos, axis=1)
        _ = sns.heatmap(scaled_genos, cmap='Blues',  # 'Greys',
                        cbar_kws = {"shrink": .5},
                        annot=cross_genos, ax=ax,
                        robust=True, square=True, mask=heatmask,
                        xticklabels=imputed_labels, yticklabels=true_labels)
        # border succesfully decoded with green
        ax.add_patch(Rectangle((3, 0), 1, 1, fill=False, edgecolor='xkcd:grass green', lw=2))
        ax.add_patch(Rectangle((4, 1), 1, 1, fill=False, edgecolor='xkcd:grass green', lw=2))
        ax.add_patch(Rectangle((5, 2), 1, 1, fill=False, edgecolor='xkcd:grass green', lw=2))
        # # border half decoded with yellow
        # ax.add_patch(Rectangle((1, 0), 1, 1, fill=False, edgecolor='xkcd:gold', lw=2))
        # ax.add_patch(Rectangle((1, 1), 1, 1, fill=False, edgecolor='xkcd:gold', lw=2))
        # ax.add_patch(Rectangle((2, 1), 1, 1, fill=False, edgecolor='xkcd:gold', lw=2))
        # ax.add_patch(Rectangle((2, 2), 1, 1, fill=False, edgecolor='xkcd:gold', lw=2))
        # # border non decoded with red
        # ax.add_patch(Rectangle((0, 0), 1, 1, fill=False, edgecolor='xkcd:rust', lw=2))
        # ax.add_patch(Rectangle((0, 1), 1, 1, fill=False, edgecolor='xkcd:rust', lw=2))
        # ax.add_patch(Rectangle((0, 2), 1, 1, fill=False, edgecolor='xkcd:rust', lw=2))
        ax.set_xlabel('Imputed GT')
        ax.set_ylabel('True GT')
        ax.set_title('Allele frequency = {}'.format(bin))
    else:
        dT = true
        dP = imputed

        confusion = np.zeros((len(imputed_genos), len(imputed_genos)))
        for y_true, y_pred in zip(dT.values, dP.values):
            cm = confusion_matrix(y_true.astype(str), y_pred.astype(str),
                                  labels=np.array(imputed_genos).astype(str))
            confusion = confusion + cm

        cross_genos = confusion[-3:]  # rectangle, no missing genotypes in the true dataset
        scaled_genos = minmax_scale(cross_genos, axis=1)
        _ = sns.heatmap(scaled_genos, cmap='Blues',  # 'Greys',
                        cbar_kws={"shrink": .5},
                        annot=cross_genos, ax=ax,
                        robust=True, square=True, mask=heatmask,
                        xticklabels=imputed_labels, yticklabels=true_labels)
        # border succesfully decoded with green
        ax.add_patch(Rectangle((3, 0), 1, 1, fill=False, edgecolor='xkcd:grass green', lw=2))
        ax.add_patch(Rectangle((4, 1), 1, 1, fill=False, edgecolor='xkcd:grass green', lw=2))
        ax.add_patch(Rectangle((5, 2), 1, 1, fill=False, edgecolor='xkcd:grass green', lw=2))
        # # border half decoded with yellow
        # ax.add_patch(Rectangle((1, 0), 1, 1, fill=False, edgecolor='xkcd:gold', lw=2))
        # ax.add_patch(Rectangle((1, 1), 1, 1, fill=False, edgecolor='xkcd:gold', lw=2))
        # ax.add_patch(Rectangle((2, 1), 1, 1, fill=False, edgecolor='xkcd:gold', lw=2))
        # ax.add_patch(Rectangle((2, 2), 1, 1, fill=False, edgecolor='xkcd:gold', lw=2))
        # # border non decoded with red
        # ax.add_patch(Rectangle((0, 0), 1, 1, fill=False, edgecolor='xkcd:rust', lw=2))
        # ax.add_patch(Rectangle((0, 1), 1, 1, fill=False, edgecolor='xkcd:rust', lw=2))
        # ax.add_patch(Rectangle((0, 2), 1, 1, fill=False, edgecolor='xkcd:rust', lw=2))
        ax.set_xlabel('Imputed GT')
        ax.set_ylabel('True GT')
        ax.set_title('Allele frequency = {}'.format('all'))
plt.savefig(os.path.join(outdir, 'genotypes_confusion.pdf'))
plt.show()