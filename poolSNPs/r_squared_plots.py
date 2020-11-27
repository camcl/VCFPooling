"""
Reference for the method:

Beagle16:
'We evaluated accuracy using the squared correlation (r2) between
the masked minor-allele dose and the imputed minor-allele dose.
The true minor-allele dose is the number of copies of the minor allele carried by an individual.
The imputed allele dose is the sum of the posterior allele probabilities for the two haplotypes of an individual.
Imputation accuracy varies with minor allele frequency, and there is little information
to estimate squared correlation for single markers when minor allele counts are low,
so we binned genotypes according to the minor allele count of the corresponding marker,
and we calculated r2 for the genotypes in each minor allele count bin.'

MaCH10:
r² between estimated and true, where true is supposed to be HWE. Since the raw values from 1KGP data comply with HWE,
they are set as the true ones.
Correlation is computed as Pearson's r² between true and imputed dosage.

Publications found do not show any p-value along with r².
"""

import os, sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr

rootdir = os.path.dirname(os.path.dirname(os.getcwd()))
sys.path.insert(0, rootdir)

from VCFPooling.poolSNPs.metrics import quality as qual
from VCFPooling.poolSNPs import dataframe as vcfdf
from VCFPooling.persotools.files import *


# Data parameters and plotting features

rW = 1000
bS = 0.02
x_data = 'maf'  # 'binned_maf_info'
x_bins = np.arange(0.0, 0.5 + bS, bS) if x_data in ['maf_info', 'maf'] \
        else np.arange(0.0, 1.0 + bS, bS)
x_bins_labels = (np.diff(x_bins) / 2) + x_bins[:-1]

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

outdir = '/home/camille/PoolImpHuman/results/20200923'  # '/home/camille/1000Genomes/src/VCFPooling/examples'
if not os.path.exists(outdir):
    os.mkdir(outdir)

true = '/home/camille/PoolImpHuman/data/20200923/IMP.chr20.snps.gt.vcf.gz'
pooled = '/home/camille/PoolImpHuman/data/20200923/IMP.chr20.pooled.snps.gt.vcf.gz'
imputed_beagle = '/home/camille/PoolImpHuman/data/20200923/IMP.chr20.pooled.imputed.beagle.vcf.gz'  # '/home/camille/PoolImpHuman/data/20200710/IMP.chr20.pooled.imputed.vcf.gz'
imputed_phaser = '/home/camille/PoolImpHuman/data/20200923/IMP.chr20.pooled.imputed.vcf.gz'

# true = '/home/camille/1000Genomes/src/VCFPooling/examples/IMP.chr20.snps.gt.vcf.gz'
# imputed_beagle = '/home/camille/1000Genomes/src/VCFPooling/examples/IMP.chr20.pooled.imputed.vcf.gz'


# Load data

dftrue = vcfdf.PandasMixedVCF(true, format='GT', indextype='chrom:pos')
maf_true = dftrue.maf
af_true = dftrue.aaf
true_bins = pd.cut(maf_true.values.flatten(), bins=x_bins, labels=x_bins_labels, include_lowest=True)
binned_maf_true = pd.DataFrame(data=true_bins, index=maf_true.index, columns=['binned_maf'])

dfpooled = vcfdf.PandasMixedVCF(pooled, format='GT', indextype='chrom:pos')
maf_pooled = dfpooled.maf
af_pooled = dfpooled.aaf
pooled_bins = pd.cut(maf_pooled.values.flatten(), bins=x_bins, labels=x_bins_labels, include_lowest=True)
binned_maf_pooled = pd.DataFrame(data=pooled_bins, index=maf_pooled.index, columns=['binned_maf'])

dfbeagle = vcfdf.PandasMixedVCF(imputed_beagle, format='GT', indextype='chrom:pos')
maf_beagle = dfbeagle.maf
af_beagle = dfbeagle.aaf

dfphaser = vcfdf.PandasMixedVCF(imputed_phaser, format='GT', indextype='chrom:pos')
maf_phaser = dfphaser.maf
af_phaser = dfphaser.aaf

dosage = {
    'true_alt': pd.Series(2 * af_true.values.flatten(), index=af_true.index, name='true_alt_dos'),
    'pooled_alt': pd.Series(2 * af_pooled.values.flatten(), index=af_pooled.index, name='pooled_alt_dos'),
    'beagle_alt': pd.Series(2 * af_beagle.values.flatten(), index=af_beagle.index, name='beagle_alt_dos'),
    'phaser_alt': pd.Series(2 * af_phaser.values.flatten(), index=af_phaser.index, name='phaser_alt_dos')
}

for k, dos in dosage.items():
    jsonf = os.path.join(outdir, '{}.json'.format(k))
    dos.to_frame().to_json(jsonf, orient='records')

dosdata = {
    'true_alt': None,
    'pooled_alt': None,
    'beagle_alt': None,
    'phaser_alt': None
}

for k in dosdata.keys():
    jsonf = os.path.join(outdir, '{}.json'.format(k))
    dosdata[k] = pd.read_json(jsonf, orient='records')

#TODO: dosage from posterior gprobs

# Correlation imputed dosage and true dosage
binned_alt_true_beagle = binned_maf_true.join([dosage['true_alt'], dosage['beagle_alt']])
corr_mat_beagle = binned_alt_true_beagle.groupby('binned_maf').corr(method='pearson')
corr_ser_beagle = corr_mat_beagle.xs('true_alt_dos', level=1)['beagle_alt_dos'].rename('r_squared')
corr_df_beagle = pd.DataFrame(data='beagle', index=corr_ser_beagle.index, columns=['dataset']).join(corr_ser_beagle)

binned_alt_true_phaser = binned_maf_true.join([dosage['true_alt'], dosage['phaser_alt']])
corr_mat_phaser = binned_alt_true_phaser.groupby('binned_maf').corr(method='pearson')
corr_ser_phaser = corr_mat_phaser.xs('true_alt_dos', level=1)['phaser_alt_dos'].rename('r_squared')
corr_df_phaser = pd.DataFrame(data='phaser', index=corr_ser_phaser.index, columns=['dataset']).join(corr_ser_phaser)

corr_df_true = pd.concat([corr_df_beagle, corr_df_phaser], axis=0).reset_index()

g = sns.lineplot(data=corr_df_true, x='binned_maf', y='r_squared', hue='dataset', lw=3)
g.set_xlabel('True minor allele frequency in study population')
g.set_ylabel('r²')
plt.title('Squared correlation between the true and the imputed allelic dosages')
plt.savefig(os.path.join(outdir, 'r_sqr_dos_true_imputed.pdf'))
plt.show()


# # Correlation imputed dosage and pooled dosage
# binned_alt_pooled_beagle = binned_maf_pooled.join([dosage['pooled_alt'], dosage['beagle_alt']])
# corr_mat_beagle = binned_alt_pooled_beagle.groupby('binned_maf').corr(method='pearson')
# corr_ser_beagle = corr_mat_beagle.xs('pooled_alt_dos', level=1)['beagle_alt_dos'].rename('r_squared')
# corr_df_beagle = pd.DataFrame(data='beagle', index=corr_ser_beagle.index, columns=['dataset']).join(corr_ser_beagle)
#
# binned_alt_pooled_phaser = binned_maf_pooled.join([dosage['pooled_alt'], dosage['phaser_alt']])
# corr_mat_phaser = binned_alt_pooled_phaser.groupby('binned_maf').corr(method='pearson')
# corr_ser_phaser = corr_mat_phaser.xs('pooled_alt_dos', level=1)['phaser_alt_dos'].rename('r_squared')
# corr_df_phaser = pd.DataFrame(data='phaser', index=corr_ser_phaser.index, columns=['dataset']).join(corr_ser_phaser)
#
# corr_df_pooled = pd.concat([corr_df_beagle, corr_df_phaser], axis=0).reset_index()
#
# g = sns.lineplot(data=corr_df_pooled, x='binned_maf', y='r_squared', hue='dataset', lw=3)
# g.set_xlabel('pooled minor allele frequency in study population')
# g.set_ylabel('r²')
# plt.title('Squared correlation between the pooled and the imputed allelic dosages')
# plt.savefig(os.path.join(outdir, 'r_sqr_dos_true_imputed.pdf'))
# plt.show()

