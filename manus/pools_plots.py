import os, sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from typing import *

rootdir = os.path.dirname(os.path.dirname(os.getcwd()))
sys.path.insert(0, rootdir)

from VCFPooling.poolSNPs.metrics import quality as qual
from VCFPooling.poolSNPs import dataframe as vcfdf
from VCFPooling.persotools.files import *

"""
RGB proportional to GL(RR|RA|AA)
"""

# TODO:
# TODO: add command-line arguments parsing
# TODO: pick relavant markers (middle-range AF and low AF) + show true-pooled-imputed --> create VCF with this marker only (calculation upspeeding)
# TODO: zoom plot for 1 block only


outdir = '/home/camille/PoolImpHuman/results/20200827'
if not os.path.exists(outdir):
    os.mkdir(outdir)

true = '/home/camille/PoolImpHuman/data/20200827/IMP.chr20.snps.gt.vcf.gz'
pooled = '/home/camille/PoolImpHuman/data/20200827/IMP.chr20.pooled.snps.gl.vcf.gz'
imputed = '/home/camille/PoolImpHuman/data/20200827/IMP.chr20.pooled.imputed.vcf.gz'


# Colouring and labelling functions

def gl_colors(gen, min_log_gl = -12):
    """

    :param gen: genotype (log-)likelihoods. Array-like.
    :return: triplet coding rgb
    """
    gen = np.array([g for g in gen])
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
    pool_labels = np.empty(var.shape[0], dtype='<U3')  # 'x/x' is 3-character long

    if af_info is None:
        af_info = np.nan
    if aaf is None:
        aaf = np.nan

    plt.rcParams['axes.titlesize'] = 12
    plots_sz = (3, 5)
    fig, axes = plt.subplots(plots_sz[0], plots_sz[1],
                             figsize=(8, 8))
    fig.suptitle('SNP 20:{0}; AF_INFO = {1:.5f}; AAF = {2:.5f}'.format(snp, af_info, aaf),
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
    # print(list(zip(pool_colors, pool_labels)))
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
    fig.suptitle('SNP 20:{0}; AF_INFO = {1:.5f}; AAF = {2:.5f}'.format(snp, af_info, aaf),
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
    # print(list(zip(pool_colors, pool_labels)))
    k = 0
    #for i in range(1, plots_sz[0] + 1):
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
    plt.show()
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
    with open(os.path.join(outdir, 'pools_patterns.snp{}.jpg'.format(snp)), "wb") as f:
        f.write(img2pdf.convert([i for i in imglist]))


# Read files

'''
20:63799     0.439497
20:68749     0.571286 --> pool #7
20:144570    0.515775
20:182013    0.520567
20:217281    0.528355
20:263660    0.519968
20:264365    0.512181  --> pools #1, #2, #4 --> row 1
20:62915126  0.006190  --> pools #3, #4, #14 relevant with ALT carrier --> row 1
'''
myvar = '20:264365'
# myvar = '20:62915126'

dftrue = vcfdf.PandasMixedVCF(true, format='GT', indextype='chrom:pos')
dfpooled = vcfdf.PandasMixedVCF(pooled, format='GL', indextype='chrom:pos')
dfimputed = vcfdf.PandasMixedVCF(imputed, format='GT', indextype='chrom:pos')

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

