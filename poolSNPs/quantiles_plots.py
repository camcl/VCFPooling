"""
Metrics across MAF range with quantiles dispersion
"""

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
bS = 0.01
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

truegt = '/home/camille/PoolImpHuman/data/20201031/IMP.chr20.snps.gt.vcf.gz'
truegl = '/home/camille/PoolImpHuman/data/20201029/IMP.chr20.snps.gl.vcf.gz'
imputed_beagle = '/home/camille/PoolImpHuman/data/20201029/IMP.chr20.pooled.imputed.vcf.gz'
imputed_phaser = '/home/camille/PoolImpHuman/data/20201031/IMP.chr20.pooled.imputed.vcf.gz'

print('\r\nData written to {}'.format(outdir))


# Function/Tools

def rollquants(dX: pd.DataFrame, dS1: pd.Series, dS2: pd.Series) -> pd.DataFrame:
    pdf1 = qual.QuantilesDataFrame(dX,
                                   dS1,
                                   bins_step=bS)
    pctY1 = pdf1.binnedX_rolling_quantilY(rollwin=rQ)
    pctY1['dataset'] = ['beagle'] * pctY1.shape[0]

    pdf2 = qual.QuantilesDataFrame(dX,
                                   dS2,
                                   bins_step=bS)
    pctY2 = pdf2.binnedX_rolling_quantilY(rollwin=rQ)
    pctY2['dataset'] = ['phaser'] * pctY2.shape[0]

    rollquants = pd.concat([pctY1, pctY2])

    return rollquants


# Load data and check

qbeaglegt = qual.QualityGT(truegt, imputed_beagle, 0, idx='chrom:pos')
qbeaglegl = qual.QualityGL(truegl, imputed_beagle, 0, idx='chrom:pos')

print('\r\n{} variants from {} samples read from {}'.format(len(qbeaglegt.trueobj.variants),
                                                            len(qbeaglegt.trueobj.samples),
                                                            os.path.basename(truegt)))
print('\r\n{} variants from {} samples read from {}'.format(len(qbeaglegt.imputedobj.variants),
                                                            len(qbeaglegt.imputedobj.samples),
                                                            os.path.basename(imputed_beagle)))
bgldiff = qbeaglegt.diff()

qphasergt = qual.QualityGT(truegt, imputed_phaser, 0, idx='chrom:pos')
qphasergl = qual.QualityGL(truegl, imputed_phaser, 0, idx='chrom:pos')

print('\r\n{} variants from {} samples read from {}'.format(len(qbeaglegl.trueobj.variants),
                                                            len(qbeaglegl.trueobj.samples),
                                                            os.path.basename(truegl)))
print('\r\n{} variants from {} samples read from {}'.format(len(qphasergl.imputedobj.variants),
                                                            len(qphasergl.imputedobj.samples),
                                                            os.path.basename(imputed_phaser)))

print('True AF INFO: \r{}'.format(qbeaglegt.trueobj.aaf))  # af_info
mafS = qbeaglegt.trueobj.maf  # maf_info

metrics = {'precision_score': {'beagle': qbeaglegt.precision,
                         'phaser': qphasergt.precision},
           'recall_score': {'beagle': qbeaglegt.recall,
                         'phaser': qphasergt.recall},
           'f1_score':  {'beagle': qbeaglegt.f1_score,
                         'phaser': qphasergt.f1_score},
           'concordance': {'beagle': qbeaglegt.concordance(),
                           'phaser': qphasergt.concordance()},
           'allelic_dos': None,
           'cross_entropy': {'beagle': qbeaglegl.cross_entropy,
                           'phaser': qphasergl.cross_entropy}
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
        yS_beagle, yS_phaser = list(d.values())
        # Compute quantiles
        print('Computing quantiles for {}'.format(metric).ljust(80, '.'))
        pctY_comp = rollquants(mafS, yS_beagle, yS_phaser)
        # Compute mean over all markers
        print('Computing means for {}'.format(metric).ljust(80, '.'))
        pctY_comp['mean'] = pctY_comp['dataset'].apply(lambda x: yS_beagle.mean() if x == 'beagle' else yS_phaser.mean())
        jsonf = dataquants[metric]
        pctY_comp.to_json(jsonf,
                          orient='records')
print(pctY_comp[pctY_comp['dataset'] == 'phaser'])


# Read processed reshaped data for plotting and draw figures

for dquant, f in dataquants.items():
    if f is not None:
        dataf = pd.read_json(f, orient='records')
        meanf = {}

        gY = sns.lineplot(data=dataf[dataf.quantiles == 0.5], x=x_data, y=dquant,
                          hue='dataset', palette="deep", linewidth=1)
        for i, dset in enumerate(['beagle', 'phaser']):
            df = dataf[dataf['dataset'] == dset]
            meanf[dset] = df['mean'].mean()
            gY.fill_between(df[df.quantiles == 1.0][x_data],
                            df[df.quantiles == 0.0][dquant],
                            df[df.quantiles == 1.0][dquant],
                            color=sns.color_palette('deep')[i],
                            alpha=0.1)
            gY.fill_between(df[df.quantiles == 0.99][x_data],
                            df[df.quantiles == 0.01][dquant],
                            df[df.quantiles == 0.99][dquant],
                            color=sns.color_palette('deep')[i],
                            alpha=0.25)
            gY.fill_between(df[df.quantiles == 0.75][x_data],
                            df[df.quantiles == 0.25][dquant],
                            df[df.quantiles == 0.75][dquant],
                            color=sns.color_palette('deep')[i],
                            alpha=0.40)
        gY.set_xlabel('True minor allele frequency in {} population'.format('study' if x_data == 'binned_maf'
                                                                            else 'main'))
        handles, labels = gY.get_legend_handles_labels()
        labels[-2] = '{} (mean = {:.5f})'.format(labels[-2], meanf['beagle'])
        labels[-1] = '{} (mean = {:.5f})'.format(labels[-1], meanf['phaser'])
        gY.legend(handles, labels)
        plt.savefig(os.path.join(outdir, '{}_percentiles_rQ={}_bS={}_xdata={}.pdf'.format(dquant, rQ, bS, x_data.lstrip('binned_'))))
        plt.show()
