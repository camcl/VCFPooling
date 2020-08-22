"""
DRAFT!!!
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

rQ = 1000
bS = 0.01
x_data = 'binned_maf'  # 'binned_maf_info'

true = '/home/camille/PoolImpHuman/data/20200722/IMP.chr20.snps.gt.vcf.gz'
imputed_beagle = '/home/camille/PoolImpHuman/data/20200722/IMP.chr20.pooled.imputed.vcf.gz'
imputed_phaser = '/home/camille/PoolImpHuman/data/20200817/IMP.chr20.pooled.imputed.vcf.gz'
# true = imputed_phaser  # test I get a concordance of 1 everywhere

qbeaglegt = qual.QualityGT(true, imputed_beagle, 0, idx='chrom:pos')

bgldiff = qbeaglegt.diff()

qphasergt = qual.QualityGT(true, imputed_phaser, 0, idx='chrom:pos')
print(qbeaglegt.trueobj.aaf)  # af_info
mafS = qbeaglegt.trueobj.maf  # maf_info

metrics = {'precision_score': {'beagle': qbeaglegt.precision,
                         'phaser': qphasergt.precision},
           # 'recall_score': {'beagle': qbeaglegt.recall,
           #               'phaser': qphasergt.recall},
           # 'f1_score':  {'beagle': qbeaglegt.f1_score,
           #               'phaser': qphasergt.f1_score},
           'concordance': {'beagle': qbeaglegt.concordance(),
                           'phaser': qphasergt.concordance()},
           'allelic_dos': None,
           'cross_entropy': None
           }

dataquants = {'precision_score': '/home/camille/PoolImpHuman/results/rolling_quantiles_precision_score.json',
              # 'recall_score': '/home/camille/PoolImpHuman/results/rolling_quantiles_recall_score.json',
              # 'f1_score': '/home/camille/PoolImpHuman/results/rolling_quantiles_f1_score.json',
              'concordance': '/home/camille/PoolImpHuman/results/rolling_quantiles_concordance.json',
              'allelic_dos': None,
              'cross_entropy': None
              }

# # roll-mean values for smoothing
# dfjoin_beagle = dX.join(dS1).sort_values(by='maf_info', axis=0)
# dfrolled_beagle = dfjoin_beagle.rolling(window=50).mean().dropna()
# dfjoin_phaser = dX.join(dS2).sort_values(by='maf_info', axis=0)
# dfrolled_phaser = dfjoin_phaser.rolling(window=50).mean().dropna()
#
# pdf_beagle = qual.QuantilesDataFrame(dfrolled_beagle[['maf_info']],
#                                      dfrolled_beagle['concordance'],
#                                      bins_step=0.005)
# pdf_phaser = qual.QuantilesDataFrame(dfrolled_phaser[['maf_info']],
#                                      dfrolled_phaser['concordance'],
#                                      bins_step=0.005)


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


for metric, d in metrics.items():
    if d is not None:
        yS_beagle, yS_phaser = list(d.values())
        pctY_comp = rollquants(mafS, yS_beagle, yS_phaser)
        jsonf = dataquants[metric]
        pctY_comp.to_json(jsonf,
                          orient='records')
print(pctY_comp[pctY_comp['dataset'] == 'phaser'])

for dquant, f in dataquants.items():
    if f is not None:
        dataf = pd.read_json(f, orient='records')
        print(dataf[dataf['dataset'] == 'phaser'])
        
        gY = sns.lineplot(data=dataf[dataf.quantiles == 0.5], x=x_data, y=dquant,
                          hue='dataset', palette="deep", linewidth=1)
        for i, dset in enumerate(['beagle', 'phaser']):
            df = dataf[dataf['dataset'] == dset]
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
        gY.set_xlabel('Estimated minor allele frequency in {} population'.format('study' if x_data == 'binned_maf'
                                                                                 else 'main'))
        
        plt.savefig('/home/camille/PoolImpHuman/results/{}_percentiles_rQ={}_bS={}_xdata={}.pdf'.format(dquant, rQ, bS, x_data.lstrip('binned_')))
        plt.show()
