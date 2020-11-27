import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, to_rgba_array, to_rgba
import sys, os
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
from collections import Counter

# force PYTHONPATH to look into the project directory for modules
rootdir = os.path.dirname(os.path.dirname(os.getcwd()))
sys.path.insert(0, rootdir)

from VCFPooling.poolSNPs import dataframe as vcfdf

"""
Plot Early Data Analysis from pooling simulation:
The baseline is OmniExpress LD STU population + HD ref PAN (missing data, genotypes distribution)
vs.
OmniExpress HD in pooled STU population + HD ref PAN
- hypergeometric distortion
- missing data rate
- genotypes het/hom_ref/hom_alt proportions
- markers counts per MAF bin

To be included in BMC Bioinformatics article

Command line usage (assuming the current directory is VCFPooling/examples)
$ python3 -u ../manus/eda.py /home/camille/PoolImpHuman/data/20200812/IMP.chr20.snps.gt.vcf.gz /home/camille/PoolImpHuman/data/20200812/IMP.chr20.missingHD.fullLD.snps.gt.vcf.gz /home/camille/PoolImpHuman/data/20200812/IMP.chr20.pooled.snps.gt.vcf.gz /home/camille/PoolImpHuman/results/20200812
"""

### COMMAND-LINE PARSING AND PARAMETERS
# parser = argparse.ArgumentParser(description='Early Data Analysis with plots'
#                                              'for the genetic structure  in a population'
#                                              'of simulated pooled genotypes')
# parser.add_argument('truegenos', metavar='trug', type=str, help='File with HD true genotypes GT', default=None)
# parser.add_argument('lowdgenos', metavar='lowdg', type=str, help='File with LD true genotypes GT and missing HD genotypes', default=None)
# parser.add_argument('pooledgenos', metavar='poog', type=str, help='File with HD pooled genotypes decoded into GT', default=None)
# parser.add_argument('outdir', metavar='outdir', type=str, help='Directory to save the plots', default=None)
# argsin = parser.parse_args()
#
# paths = {'gt': {
#     'true': argsin.truegenos,
#     'lowd': argsin.lowdgenos,
#     'pooled': argsin.pooledgenos}
# }

paths = {'gt': {
    'true': '/home/camille/PoolImpHuman/data/20200812/IMP.chr20.snps.gt.vcf.gz',
    'lowd': '/home/camille/PoolImpHuman/data/20200812/IMP.chr20.missingHD.fullLD.snps.gt.vcf.gz',
    'pooled': '/home/camille/PoolImpHuman/data/20200812/IMP.chr20.pooled.snps.gt.vcf.gz'}
}

outdir = '/home/camille/PoolImpHuman/results/20200812'  # argsin.outdir
if not os.path.exists(outdir):
    os.mkdir(outdir)
print('\r\nFigures and tables will be saved in {}'.format(outdir).ljust(80, '.'))


# Read and process data

dftrue = vcfdf.PandasMixedVCF(paths['gt']['true'], format='GT')
dflowd = vcfdf.PandasMixedVCF(paths['gt']['lowd'], format='GT')
dfpool = vcfdf.PandasMixedVCF(paths['gt']['pooled'], format='GT')

af_data = 'aaf'
df0 = dftrue.concatcols([dftrue.af_info, dftrue.missing_rate, dftrue.aaf, dftrue.het_rate, dftrue.hom_alt_rate, dftrue.hom_ref_rate])
df1 = dflowd.concatcols([dftrue.af_info, dftrue.aaf.rename(columns={'aaf': 'true_af'}), dflowd.missing_rate, dflowd.aaf])
df2 = dfpool.concatcols([dftrue.af_info, dftrue.aaf.rename(columns={'aaf': 'true_af'}), dfpool.missing_rate, dfpool.aaf])
df3 = dflowd.concatcols([dflowd.af_info, dftrue.aaf.rename(columns={'aaf': 'true_af'}), dflowd.het_rate, dflowd.hom_alt_rate, dflowd.hom_ref_rate])
df4 = dfpool.concatcols([dfpool.af_info, dftrue.aaf.rename(columns={'aaf': 'true_af'}), dfpool.het_rate, dfpool.hom_alt_rate, dfpool.hom_ref_rate])


print(df1.head(10))
print(df2.head(10))

print('\r\nOverall missing rates:')
totldmiss = dflowd.missing_rate.mean()
totpoolmiss = dfpool.missing_rate.mean()
print('From LD to HD: {}'.format(totldmiss))
print('In pooled data: {}'.format(totpoolmiss))


# Plots' layout params

pts_sz = 10  # scatter dots size
pts_alpha = 1.0  # transparency
figsize = 8
labfont = 10

true_genos = [0.0, 1.0, 2.0]
true_labels = ['0/0', '0/1', '1/1']
pooled_genos = [0.0, 1.0, 2.0, -0.5, 0.5, -1.0]
pooled_labels = ['0/0', '0/1', '1/1', '0/.', './1', './.']
genocolors = ['#047495', '#00035b', '#748b97',  # full GT
              '#dbb40c', '#c65102', '#80013f'  # missing GT
              ]
genotricmap = ListedColormap([to_rgba(co) for co in genocolors[:3]], name='geno_tri_cmap')

x_bins = [0.0, 0.02, 0.04, 0.06, 0.1, 0.2, 0.4, 0.6, 0.8, 0.9, 0.94, 0.96, 0.98, 1.0]
lab_bins = [0.01, 0.03, 0.05, 0.08, 0.15, 0.3, 0.5, 0.7, 0.85, 0.92, 0.95, 0.97, 0.99]


# Basic statistics for LD and HD data sets per bin

# Markers counts

cnts = dict([(lab, 0) for lab in lab_bins])
ldbin = pd.cut(dflowd.aaf.values.squeeze(), x_bins, labels=lab_bins, include_lowest=True)
cnts.update(Counter(ldbin.dropna()))
ldcnts = cnts

cnts = dict([(lab, 0) for lab in lab_bins])
hdbin = pd.cut(dftrue.aaf.values.squeeze(), x_bins, labels=lab_bins, include_lowest=True)
cnts.update(Counter(hdbin.dropna()))
hdcnts = cnts

dfbincnts = pd.DataFrame.from_records([ldcnts, hdcnts])
dfbincnts.index = ['LD', 'HD']
dfbincnts['Total'] = dfbincnts.sum(axis=1)
dfbincnts.to_latex(buf=os.path.join(outdir, 'markers-LD-HD-bin-counts.tex'),
                   sparsify=True,
                   multirow=True,
                   caption='SNPs counts per AAF bin: Columns labels are the central values of each interval',
                   label='tab:markers-LD-HD-bin-counts')

# Counts for rare variants only
rarecnts = dfbincnts[dfbincnts.columns[:3]].sum(axis=1) + dfbincnts[dfbincnts.columns[-4:-1]].sum(axis=1)
# LD    1972
# HD    20833


# Missing counts

ldbin = pd.cut(df1.true_af.values.squeeze(), x_bins, labels=lab_bins, include_lowest=True)
df1['AF-bin'] = ldbin
ldhdmiss = df1.groupby(['AF-bin'])['missing_rate'].mean().to_frame()

pooledbin = pd.cut(df2.true_af.values.squeeze(), x_bins, labels=lab_bins, include_lowest=True)
df2['AF-bin'] = pooledbin
pooledmiss = df2.groupby(['AF-bin'])['missing_rate'].mean().to_frame()

binmiss = ldhdmiss.join(pooledmiss, how='right', lsuffix='_LDHD', rsuffix='_pooledHD').transpose()
binmiss.columns = binmiss.columns.to_list()
binmiss['Total'] = np.concatenate([totldmiss.values, totpoolmiss.values])
binmiss.index = ['LD + HD', 'pooled HD']
binmiss.to_latex(buf=os.path.join(outdir, 'missing-LDHD-pooled-bin.tex'),
                 sparsify=True,
                 multirow=True,
                 float_format="%.3f",
                 caption='Proportion of missing genotypes per AAF bin: Columns labels are the central values of each interval',
                 label='tab:missing-LDHD-pooled-bin')
