"""
Cast the results into LaTeX-formatted tables instead of plots"""

import os, sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from typing import *

rootdir = os.path.dirname(os.path.dirname(os.getcwd()))
sys.path.insert(0, rootdir)

from VCFPooling.poolSNPs.metrics import quality as qual
from VCFPooling.persotools.files import *


# Data parameters and plotting features

rQ = 1000
bS = 0.1
x_data = 'binned_maf'  # 'binned_maf_info'

sns.set(rc={'figure.figsize': (10, 8)})
sns.set_style('whitegrid')
dash_styles = [
    (1, 1),
    (3, 1, 1.5, 1),
    (5, 1, 1, 1),
    (5, 1, 2, 1, 2, 1),
    (2, 2, 3, 1.5),
    (1, 2.5, 3, 1.2),
    "",
    (4, 1.5),
]


# Configure data/plots paths

outdir = '/home/camille/PoolImpHuman/results/20201031'
if not os.path.exists(outdir):
    os.mkdir(outdir)

truegt = '/home/camille/PoolImpHuman/data/20200827/sHG01063.IMP.chr20.snps.gt.vcf.gz'
truegl = '/home/camille/PoolImpHuman/data/20200827/sHG01063.IMP.chr20.snps.gl.vcf.gz'
imputed_phaser1 = '/home/camille/PoolImpHuman/data/20201028/sHG01063.IMP.chr20.pooled.imputed.vcf.gz'
imputed_phaser2 = '/home/camille/PoolImpHuman/data/20201031/sHG01063.IMP.chr20.snps.gl.full.postgenos.vcf.gz'


# Function/Tools

def rollquants(dX: pd.DataFrame, dS1: pd.Series, dS2: pd.Series) -> pd.DataFrame:
    pdf1 = qual.QuantilesDataFrame(dX,
                                   dS1,
                                   bins_step=bS)
    pctY1 = pdf1.binnedX_rolling_quantilY(rollwin=rQ)
    pctY1['dataset'] = ['phaser1'] * pctY1.shape[0]

    pdf2 = qual.QuantilesDataFrame(dX,
                                   dS2,
                                   bins_step=bS)
    pctY2 = pdf2.binnedX_rolling_quantilY(rollwin=rQ)
    pctY2['dataset'] = ['phaser2'] * pctY2.shape[0]

    rollquants = pd.concat([pctY1, pctY2])

    return rollquants


# Load data

qphaser1gt = qual.QualityGT(truegt, imputed_phaser1, 0, idx='chrom:pos')
qphaser1gl = qual.QualityGL(truegl, imputed_phaser1, 0, idx='chrom:pos')

bgldiff = qphaser1gt.diff()

qphaser2gt = qual.QualityGT(truegt, imputed_phaser2, 0, idx='chrom:pos')
qphaser2gl = qual.QualityGL(truegl, imputed_phaser2, 0, idx='chrom:pos')

print(qphaser1gt.trueobj.aaf)  # af_info
mafS = qphaser1gt.trueobj.maf  # maf_info

metrics = {'precision_score': {'phaser1': qphaser1gt.precision,
                         'phaser2': qphaser2gt.precision},
           'recall_score': {'phaser1': qphaser1gt.recall,
                         'phaser2': qphaser2gt.recall},
           'f1_score':  {'phaser1': qphaser1gt.f1_score,
                         'phaser2': qphaser2gt.f1_score},
           'concordance': {'phaser1': qphaser1gt.concordance(),
                           'phaser2': qphaser2gt.concordance()},
           'allelic_dos': None,
           'cross_entropy': {'phaser1': qphaser1gl.cross_entropy,
                           'phaser2': qphaser2gl.cross_entropy}
           }

dataquants = {'precision_score': os.path.join(outdir, 'rolling_quantiles_precision_score.json'),
              'recall_score': os.path.join(outdir, 'rolling_quantiles_recall_score.json'),
              'f1_score': os.path.join(outdir, 'rolling_quantiles_f1_score.json'),
              'concordance': os.path.join(outdir, 'rolling_quantiles_concordance.json'),
              'allelic_dos': None,
              'cross_entropy': os.path.join(outdir, 'rolling_quantiles_cross_entropy.json')
              }


# Process and write data

for metric, d in metrics.items():
    if d is not None:
        yS_phaser1, yS_phaser2 = list(d.values())
        # Compute quantiles
        print('Computing quantiles for {}'.format(metric).ljust(80, '.'))
        pctY_comp = rollquants(mafS, yS_phaser1, yS_phaser2)
        # Compute mean over all markers
        print('Computing means for {}'.format(metric).ljust(80, '.'))
        pctY_comp['mean'] = pctY_comp['dataset'].apply(lambda x: yS_phaser1.mean() if x == 'phaser1' else yS_phaser2.mean())
        jsonf = dataquants[metric]
        pctY_comp.to_json(jsonf,
                          orient='records')
print(pctY_comp[pctY_comp['dataset'] == 'phaser2'])


# Read processed reshaped data for plotting and draw figures

qframes = {}
for dquant, f in dataquants.items():
    if f is not None:
        dataf = pd.read_json(f, orient='records')
        meanf = {}
        qframes[dquant] = dataf

qfconc = qframes['concordance']
# qfconc[qfconc['dataset'] == 'phaser1'][qfconc['quantiles'] == 0.50]
qfconc.set_index(['binned_maf', 'quantiles'], inplace=True)
qfconc.index.names=['MAF bin', 'quantile']
qfconc1 = qfconc[qfconc['dataset'] == 'phaser1']
qfconc1['concordance'].to_latex(buf=os.path.join(outdir, 'rolling_quantiles_concordance.tex'),
                                sparsify=True,
                                multirow=True,
                                caption='Phaser not logged GL')
