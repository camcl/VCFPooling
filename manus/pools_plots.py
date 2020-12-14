#!/usr/bin/python3
# -*- coding: utf-8 -*-

import os, sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import pysam

from typing import *

rootdir = os.path.dirname(os.path.dirname(os.getcwd()))
sys.path.insert(0, rootdir)

from VCFPooling.poolSNPs import pybcf
from VCFPooling.poolSNPs import dataframe as vcfdf
from VCFPooling.persotools.files import *

"""
For a given variant, plot the pooling blocks as coloured square matrices with 
genotypes displayed each with RGB color proportional to GL(RR|RA|AA).
Useful to get a picture of how genotypes are processed from the true ones into pooled and finally imputed ones

Usage example: (command-line parsing not working yet)
VCFPooling/manus$ python3 -u /home/camille/PoolImpHuman/data/20200827/IMP.chr20.snps.gt.vcf.gz /home/camille/PoolImpHuman/data/20200827/IMP.chr20.pooled.snps.gl.vcf.gz /home/camille/PoolImpHuman/data/20200827/IMP.chr20.pooled.imputed.vcf.gz /home/camille/PoolImpHuman/results/20200827 20 62915126   
"""

# TODO: docstrings and annotations

# # Parse command-line arguments
#
# parser = argparse.ArgumentParser(description="Pooling blocks as coloured square matrices with RGB colors for samples' genotypes")
# parser.add_argument('truefile', metavar='truef', type=str, help='File with true genotypes GT', default=None)
# parser.add_argument('pooledfile', metavar='imputedf', type=str, help='File with pooled genotypes decoded into GL', default=None)
# parser.add_argument('imputedfile', metavar='imputedf', type=str, help='File with imputed genotypes decoded into GT', default=None)
# parser.add_argument('outdir', metavar='outdir', type=str, help='Directory to save the plots', default=None)
# parser.add_argument('chromosome', metavar='chrom', type=str, help='Chromosome number', default='20')
# parser.add_argument('position', metavar='pos', type=str, help='Target variant position', default=None)
#
# argsin = parser.parse_args()
#
# truef = argsin.truefile
# pooledf = argsin.pooledfile
# imputedf = argsin.imputedfile
# outdir = argsin.outdir
# chr = argsin.chromosome
# pos = argsin.position

'''
20:264365    0.512181  --> pools #1, #2, #4 --> row 1
20:62915126  0.006190  --> pools #3, #4, #14 relevant with ALT carrier --> row 1
'''

chr, pos = '20', '62915126'
myvar = '{}:{}'.format(chr, pos)


# Paths to files to read genotypes from

outdir = '/home/camille/PoolImpHuman/results/20200827'
if not os.path.exists(outdir):
    os.mkdir(outdir)

truef = '/home/camille/PoolImpHuman/data/20200827/IMP.chr20.snps.gt.vcf.gz'
pooledf = '/home/camille/PoolImpHuman/data/20200827/IMP.chr20.pooled.snps.gl.vcf.gz'
imputedf = '/home/camille/PoolImpHuman/data/20200827/IMP.chr20.pooled.imputed.vcf.gz'

for f in [truef, pooledf, imputedf]:
    pybcf.index(f, os.path.dirname(f))


# Create files for the chosen SNP only (speed up processing then

truevar = pybcf.view_one_variant(truef, pos, chr, wd=outdir)
pooledvar = pybcf.view_one_variant(pooledf, pos, chr, wd=outdir)
imputedvar = pybcf.view_one_variant(imputedf, pos, chr, wd=outdir)


# Colouring and labelling functions

def gl_colors(gen: np.ndarray, min_log_gl: float = -12.) -> np.ndarray:
    """
    Convert log-GL values to not logged ones, used as a RGB color triplet
    :param gen: genotype (log-)likelihoods. Array-like.
    :param min_log_gl: float (negative) indicating the minimum value used instead of log(GP=0.0)
    :return: triplet coding RGB channels
    """
    gen = np.array([g for g in gen])  # forces GL-tuple to GL-array conversion
    tenpow = np.vectorize(lambda x: 0.0 if x == min_log_gl else pow(10, x))
    rgb = tenpow(gen)
    return rgb


def map_gt_gl(gt_val):
    # gt_val from trinary encoding, true or imputed VCF file
    if gt_val == 2.0:
        gl_val = np.array([0., 0., 1.])
    elif gt_val == 1.0:
        gl_val = np.array([0., 1., 0.])
    elif gt_val == 0.0:
        gl_val = np.array([1., 0., 0.])
    else:  # missing data
        gl_val = np.array([1/3, 1/3, 1/3])

    return gl_val


def gt_colors(gen):
    """

    :param gen: raw genotypes with phase. Array-like.
    :return: triplet coding rgb
    """

    rgb = np.apply_along_axis(map_gt_gl, axis=-1, arr=gen)
    return rgb


def gl_labels(gen):
    """
    Returns the most likely genotype (slots labels of the pooling blocks).
    :param gen: triplet for genotype likelihoods
    :return:
    """
    if gen[0] == 1.0:
        gt = '0/0'
    elif gen[2] == 1.0:
        gt = '1/1'
    elif gen[1] == 1.0:
        gt = '1/0'
    elif round(gen[1], 1) == 0.5 and round(gen[2], 1) == 0.5:
        gt = '1/.'
    elif round(gen[0], 1) == 0.5 and round(gen[1], 1) == 0.5:
        gt = '0/.'
    else:
        gt = './.'
    return gt


def gt_labels(gen):
    """
    Returns the label genotype (slots labels of the pooling blocks).
    :param gen: hexa or trinary encoding
    :return:
    """
    if gen == 0.0:
        gt = '0/0'
    elif gen == 2.0:
        gt = '1/1'
    elif gen == 1.0:
        gt = '1/0'
    elif gen == 0.5:
        gt = '1/.'
    elif gen == -0.5:
        gt = '0/.'
    else:
        gt = './.'
    return gt


def plot_pools(var, gtgl, step, snp, af_info=None, aaf=None):
    """

    :param var: trinary encoded genotypes or log GL
    :return:
    """
    pool_rgb = np.zeros((var.shape[0], 3), dtype=float)  # 3 color channels
    pool_labels = np.empty(var.shape[0], dtype='<U3')  # GT 'x/x' is 3-character long

    if af_info is None:
        af_info = np.nan
    if aaf is None:
        aaf = np.nan

    plt.rcParams['axes.titlesize'] = 12
    plots_sz = (3, 5)
    fig, axes = plt.subplots(plots_sz[0], plots_sz[1],
                             figsize=(8, 8))
    fig.suptitle('{0} data: SNP {1}; AF_INFO = {2:.5f}; AAF = {3:.5f}'.format(step.capitalize(), snp, af_info, aaf),
                 fontsize=12)

    if gtgl.upper() == 'GL':
        pool_rgb = np.apply_along_axis(gl_colors, -1, var)  # gl_colors(var[:, np.newaxis])
        for idx, g in enumerate(pool_rgb):  # unlogged GL colors
            pool_labels[idx] = gl_labels(g)
        pool_rgb = pool_rgb.reshape(15, 16, 3)  # block-wise reshaping
        pool_labels = pool_labels.reshape(15, 16)
    if gtgl.upper() == 'GT':
        pool_rgb = gt_colors(var[:, np.newaxis])
        for idx, g in enumerate(var):
            pool_labels[idx] = gt_labels(g)
        pool_rgb = pool_rgb.reshape(15, 16, 3)  # block-wise reshaping
        pool_labels = pool_labels.reshape(15, 16)

    pool_colors = np.multiply(pool_rgb, np.broadcast_to([0.5, 0.8, 0.5], pool_rgb.shape))  # modify RGB colors

    k = 0
    for i in range(1, plots_sz[0] + 1):
        for j in range(1, plots_sz[1] + 1):
            k += 1
            axes[i - 1, j - 1].imshow(pool_colors[k-1, :, :].reshape((4, 4, 3)),
                                      cmap='plasma')
            axes[i - 1, j - 1].set_title('Pool #{}'.format(k), pad=2)
            axes[i - 1, j - 1].set_xticks(np.arange(4) + 0.5, minor=True)
            axes[i - 1, j - 1].set_yticks(np.arange(4) + 0.5, minor=True)
            axes[i - 1, j - 1].tick_params(which="both",
                                           axis='both',
                                           bottom=False,
                                           left=False,
                                           labelsize=8,
                                           length=1,
                                           pad=0.3)
            axes[i - 1, j - 1].grid(which='minor', axis='both',
                                    color="w", linestyle='-', linewidth=0.25)
            # remnove borders
            axes[i - 1, j - 1].spines['top'].set_visible(False)
            axes[i - 1, j - 1].spines['right'].set_visible(False)
            axes[i - 1, j - 1].spines['bottom'].set_visible(False)
            axes[i - 1, j - 1].spines['left'].set_visible(False)
            tx_i_j = pool_labels[k - 1, :].reshape((4, 4))
            for m in range(4):
                for n in range(4):
                    axes[i - 1, j - 1].text(n, m, tx_i_j[m, n], ha="center", va="center", color="w", fontsize=10)
    plt.autoscale()
    fig.tight_layout()
    plt.show()
    plt.savefig(os.path.join(outdir, 'pools_patterns.{}.snp{}.jpg'.format(step, snp)),
                dpi=500)


def plot_5_first_pools(var, gtgl, step, snp, af_info=None, aaf=None):
    """
    First 5 pools only
    :param var: trinary encoded genotypes or log GL
    :return:
    """
    pool_rgb = np.zeros((var.shape[0], 3), dtype=float)  # 3 color channels
    pool_labels = np.empty(var.shape[0], dtype='<U3')  # 'x/x' is 3-character long

    if af_info is None:
        af_info = np.nan
    if aaf is None:
        aaf = np.nan

    plt.rcParams['axes.titlesize'] = 12
    plots_sz = (1, 5)
    fig, axes = plt.subplots(plots_sz[0], plots_sz[1],
                             figsize=(8, 2.5))
    fig.suptitle('{0} data: SNP {1}; AF_INFO = {2:.5f}; AAF = {3:.5f}'.format(step.capitalize(), snp, af_info, aaf),
                 fontsize=12)

    if gtgl.upper() == 'GL':
        pool_rgb = np.apply_along_axis(gl_colors, -1, var)  # gl_colors(var[:, np.newaxis])
        for idx, g in enumerate(pool_rgb):  # unlogged GL colors
            pool_labels[idx] = gl_labels(g)
        pool_rgb = pool_rgb.reshape(15, 16, 3)  # block-wise reshaping
        pool_labels = pool_labels.reshape(15, 16)
    if gtgl.upper() == 'GT':
        pool_rgb = gt_colors(var.values[:, np.newaxis])
        for idx, g in enumerate(var):
            pool_labels[idx] = gt_labels(g)
        pool_rgb = pool_rgb.reshape(15, 16, 3)  # block-wise reshaping
        pool_labels = pool_labels.reshape(15, 16)

    pool_colors = np.multiply(pool_rgb, np.broadcast_to([0.5, 0.8, 0.5], pool_rgb.shape))  # modify RGB colors

    k = 0
    for j in range(1, plots_sz[1] + 1):
        k += 1
        axes[j - 1].imshow(pool_colors[k-1, :, :].reshape((4, 4, 3)),
                           cmap='plasma')
        axes[j - 1].set_title('Pool #{}'.format(k), pad=2)
        axes[j - 1].set_xticks(np.arange(4) + 0.5, minor=True)
        axes[j - 1].set_yticks(np.arange(4) + 0.5, minor=True)
        axes[j - 1].tick_params(which="both",
                                axis='both',
                                bottom=False,
                                left=False,
                                labelsize=8,
                                length=1,
                                pad=0.3)
        axes[j - 1].grid(which='minor', axis='both',
                         color="w", linestyle='-', linewidth=0.25)
        # remnove borders
        axes[j - 1].spines['top'].set_visible(False)
        axes[j - 1].spines['right'].set_visible(False)
        axes[j - 1].spines['bottom'].set_visible(False)
        axes[j - 1].spines['left'].set_visible(False)
        tx_i_j = pool_labels[k - 1, :].reshape((4, 4))
        for m in range(4):
            for n in range(4):
                axes[j - 1].text(n, m, tx_i_j[m, n], ha="center", va="center", color="w", fontsize=10)
    plt.autoscale()
    fig.tight_layout()
    # plt.show()
    plt.savefig(os.path.join(outdir, 'pools_patterns.{}.snp{}.jpg'.format(step, snp)),
                dpi=500)


def before_after_pooling(snp):
    """
    Combine pools overview in a PDF document
    :param snp:
    :param sz:
    :return:
    """
    import img2pdf
    imglist = [os.path.join(outdir, 'pools_patterns.true.snp{}.jpg'.format(snp)),
               os.path.join(outdir, 'pools_patterns.pooled.snp{}.jpg'.format(snp)),
               os.path.join(outdir, 'pools_patterns.imputed.snp{}.jpg'.format(snp))]
    with open(os.path.join(outdir, 'pools_patterns.snp{}.pdf'.format(snp)), "wb") as f:
        f.write(img2pdf.convert([i for i in imglist]))


# Read files, extract relevant informations for plotting and plot pooling blocks

dftrue = vcfdf.PandasMixedVCF(truevar, format='GT', indextype='chrom:pos')
dfpooled = vcfdf.PandasMixedVCF(pooledvar, format='GL', indextype='chrom:pos')
dfimputed = vcfdf.PandasMixedVCF(imputedvar, format='GT', indextype='chrom:pos')

genos_true = dftrue.trinary_encoding().loc[myvar]
af_info = float(dftrue.af_info.loc[myvar].values)
aaf_true = float(dftrue.aaf.loc[myvar].values)

genos_pooled = dfpooled.genotypes().loc[myvar].values
aaf_pooled = float(dfpooled.aaf.loc[myvar].values)

genos_imputed = dfimputed.trinary_encoding().loc[myvar]
aaf_imputed = float(dfimputed.aaf.loc[myvar].values)

# plot_pools(genos_true, 'GT', 'true', myvar, af_info=af_info, aaf=aaf_true)  # OK
# plot_pools(genos_pooled, 'GL', 'pooled', myvar, af_info=af_info, aaf=None)  # OK
# plot_pools(genos_imputed, 'GT', 'imputed', myvar, af_info=af_info, aaf=aaf_imputed)  # OK

plot_5_first_pools(genos_true, 'GT', 'true', myvar, af_info=af_info, aaf=aaf_true)  # OK
plot_5_first_pools(genos_pooled, 'GL', 'pooled', myvar, af_info=af_info, aaf=None)  # OK
plot_5_first_pools(genos_imputed, 'GT', 'imputed', myvar, af_info=af_info, aaf=aaf_imputed)  # OK
before_after_pooling(myvar)
