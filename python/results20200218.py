import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.stats import *
import numpy as np
import pandas as pd
from cyvcf2 import VCF

from scripts.VCFPooling.poolSNPs.alleles import alleles_tools as alltls
from scripts.VCFPooling.poolSNPs.metrics import quality
from scripts.VCFPooling.poolSNPs import parameters as prm

from persotools.files import *

"""
Plot pooling evaluation:
- hypergeometric distortion
- missing data rate
- genotypes het/hom_ref/hom_alt proportions

To be included in BMC Bioinformatics article
"""

# Example with results from Beagle on 10000 markers: 240 samples imputed in 1 batch or by-by-one and merged
paths = {'gt': {
    'true': '/home/camille/1000Genomes/data/gt/stratified/IMP.chr20.snps.gt.chunk10000.vcf.gz',
    'pooled': '/home/camille/1000Genomes/data/gt/stratified//IMP.chr20.pooled.snps.gt.chunk10000.vcf.gz'}
}

dftrue = alltls.PandasVCF(paths['gt']['true'])
dfpool = alltls.PandasVCF(paths['gt']['pooled'])

df1 = dftrue.concatcols([dftrue.af_info, dftrue.missing_rate, dftrue.aaf])
df2 = dfpool.concatcols([dfpool.af_info, dfpool.missing_rate, dfpool.aaf])
df3 = dftrue.concatcols([dftrue.af_info, dftrue.het_rate, dftrue.hom_alt_rate, dftrue.hom_ref_rate])
df4 = dfpool.concatcols([dfpool.af_info, dfpool.het_rate, dfpool.hom_alt_rate, dfpool.hom_ref_rate])

print(df1.head())
plt.rcParams["figure.figsize"] = [6*2, 5*2]
fig, axes = plt.subplots(2, 2)

df1.plot.scatter('af_info', 'missing_rate', ax=axes[0, 0], s=0.7, c='tab:orange', label='missing_rate')
axes[0, 0].set_ylim(-0.01, 1.0)
ax002 = axes[0, 0].twinx()
ax002.set_ylim(-0.01, 1.0)
df1.plot.scatter('af_info', 'aaf', ax=ax002, s=0.7, c='tab:blue', label='aaf ')
axes[0, 0].set_title('true data')

df2.plot.scatter('af_info', 'missing_rate', ax=axes[0, 1], s=0.7, c='tab:orange', label='missing_rate')
axes[0, 1].set_ylim(-0.01, 1.0)
ax012 = axes[0, 1].twinx()
ax012.set_ylim(-0.01, 1.0)
df2.plot.scatter('af_info', 'aaf', ax=ax012, s=0.7, c='tab:blue', label='aaf ')
axes[0, 1].set_title('pooled data')

df3.plot.scatter('af_info', 'het_rate', ax=axes[1, 0], s=0.7, c='tab:olive', label='heterozygous')
df3.plot.scatter('af_info', 'hom_alt_rate', ax=axes[1, 0], s=0.7, c='tab:green', label='homozygous_alt')
df3.plot.scatter('af_info', 'hom_ref_rate', ax=axes[1, 0], s=0.7, c='tab:brown', label='homozygous_ref')
axes[1, 0].set_ylabel('genotype proportion')
axes[1, 0].set_title('true data')

df4.plot.scatter('af_info', 'het_rate', ax=axes[1, 1], s=0.7, c='tab:olive', label='heterozygous')
df4.plot.scatter('af_info', 'hom_alt_rate', ax=axes[1, 1], s=0.7, c='tab:green', label='homozygous_alt')
df4.plot.scatter('af_info', 'hom_ref_rate', ax=axes[1, 1], s=0.7, c='tab:brown', label='homozygous_ref')
axes[1, 1].set_ylabel('genotype proportion')
axes[1, 1].set_title('pooled data')

plt.show()
