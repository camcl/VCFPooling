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
Plot imputation performance in:
- precision
- accuracy
- recall,
- f1_score
- allele dosage
- Pearson's correlation
- cross-entropy
"""

convert_files = False

# Example with results from Phaser and Beagle on 10000 markers
paths = {'phaser': {
    'true': '/home/camille/1000Genomes/data/gt/stratified/IMP.chr20.snps.gt.chunk10000.vcf.gz',
    'imputed': '/home/camille/1000Genomes/data/gl/gl_adaptive/phaser/IMP.chr20.pooled.imputed.gt.chunk10000.vcf.gz'},
    'phasergl': {
    'true': '/home/camille/1000Genomes/data/gl/IMP.chr20.snps.gl.chunk10000.vcf',
    'imputed': '/home/camille/1000Genomes/data/gl/gl_adaptive/phaser/IMP.chr20.pooled.imputed.gl.chunk10000.vcf'},
    'beagle': {
        'true': '/home/camille/1000Genomes/data/gt/stratified/IMP.chr20.snps.gt.chunk10000.vcf.gz',
        'imputed': '/home/camille/1000Genomes/data/gl/gl_adaptive/chunk10000_20190725/IMP.chr20.pooled.imputed.gt.chunk10000.vcf.gz'},
    'beaglegl': {
        'true': '/home/camille/1000Genomes/data/gl/IMP.chr20.snps.gl.chunk10000.vcf',
        'imputed': '/home/camille/1000Genomes/data/gl/gl_adaptive/chunk10000_20190725/IMP.chr20.pooled.beagle2.gl.chunk10000.corr.vcf.gz'}
}

# Convert GT files to GL (Phaser)
if convert_files:
    alltls.file_likelihood_converter(paths['beagle']['true'], paths['beaglegl']['true'], func=alltls.bin_gl_converter)
    # easier for cross entropies not to log gl
    alltls.file_likelihood_converter(paths['phaser']['imputed'], paths['phasergl']['imputed'], func=alltls.bin_gl_converter)

qbeaglegl = quality.QualityGL(paths['beaglegl']['true'], paths['beaglegl']['imputed'], 0, idx='id')
messbeagle = qbeaglegl.cross_entropy

qphasergl = quality.QualityGL(paths['phasergl']['true'], paths['phasergl']['imputed'], 0, fmt='GL', idx='chrom:pos')
messphaser = qphasergl.cross_entropy

# dfaf = alltls.PandasVCF(paths['beagle']['true'])
# dfmess = mess.to_frame()
# dfmess = dfmess.join(dfaf.af_info)
# dfmess.plot.scatter('af_info', 'cross_entropy')
# plt.show()

qphaser = quality.QualityGT(*paths['phaser'].values(), 0, idx='chrom:pos')
tabphaser = pd.concat([qphaser.concordance(),
                       qphaser.trueobj.af_info,
                       qphaser.pearsoncorrelation(),
                       qphaser.precision,
                       qphaser.accuracy,
                       qphaser.recall,
                       qphaser.f1_score], axis=1)
tabphaser = tabphaser.join(messphaser)
dosphaser = qphaser.alleledosage()

qbeagle = quality.QualityGT(*paths['beagle'].values(), 0, idx='id')
tabbeagle = pd.concat([qbeagle.concordance(),
                       qbeagle.trueobj.af_info,
                       qbeagle.pearsoncorrelation(),
                       qbeagle.precision,
                       qbeagle.accuracy,
                       qbeagle.recall,
                       qbeagle.f1_score], axis=1)
tabbeagle = tabbeagle.join(messbeagle)
dosbeagle = qbeagle.alleledosage()

print(tabphaser.head())
print(tabbeagle.head())

plt.rcParams["figure.figsize"] = [5*2, 4*7]
fig, axes = plt.subplots(7, 2)

tabbeagle.plot.scatter('af_info', 'precision_score', ax=axes[0, 0], s=0.7, label='beagle')
axes[0, 0].set_ylim(0.0, 1.0)
tabbeagle.plot.scatter('af_info', 'accuracy_score', ax=axes[1, 0], s=0.7, label='beagle')
axes[1, 0].set_ylim(0.0, 1.0)
tabbeagle.plot.scatter('af_info', 'concordance', ax=axes[2, 0], s=0.7, label='beagle')
axes[2, 0].set_ylim(0.0, 1.2)
tabbeagle.plot.scatter('af_info', 'f1_score', ax=axes[3, 0], s=0.7, label='beagle')
axes[3, 0].set_ylim(0.0, 1.0)
tabbeagle.plot.scatter('af_info', 'r_squared', ax=axes[4, 0], s=0.7, label='beagle')
axes[4, 0].set_ylim(-0.4, 1.0)
tabbeagle.plot.scatter('af_info', 'cross_entropy', ax=axes[5, 0], s=0.7, label='beagle')
#axes[5, 0].set_ylim(0.0, 1.1)
axes[6, 0].scatter(dosbeagle[0], dosbeagle[1], s=0.7, label='beagle')
axes[6, 0].set_xlabel('true allele dosage')
axes[6, 0].set_ylabel('imputed allele dosage')
axes[6, 0].set_ylim(0.0, 2.0)

tabphaser.plot.scatter('af_info', 'precision_score', ax=axes[0, 1], s=0.7, label='phaser')
axes[0, 1].set_ylim(0.0, 1.0)
tabphaser.plot.scatter('af_info', 'accuracy_score', ax=axes[1, 1], s=0.7, label='phaser')
axes[1, 1].set_ylim(0.0, 1.0)
tabphaser.plot.scatter('af_info', 'concordance', ax=axes[2, 1], s=0.7, label='phaser')
axes[2, 1].set_ylim(0.0, 1.2)
tabphaser.plot.scatter('af_info', 'f1_score', ax=axes[3, 1], s=0.7, label='phaser')
axes[3, 1].set_ylim(0.0, 1.0)
tabphaser.plot.scatter('af_info', 'r_squared', ax=axes[4, 1], s=0.7, label='phaser')
axes[4, 1].set_ylim(-0.4, 1.0)
tabphaser.plot.scatter('af_info', 'cross_entropy', ax=axes[5, 1], s=0.7, label='phaser')
#axes[5, 1].set_ylim(0.0, 1.1)
axes[6, 1].scatter(dosphaser[0], dosphaser[1], s=0.7, label='phaser')
axes[6, 1].set_xlabel('true allele dosage')
axes[6, 1].set_ylabel('imputed allele dosage')
axes[6, 1].set_ylim(0.0, 2.0)
for ax in axes.flatten()[:-2]:
    ax.set_xlabel('true alternate allele frequency')
plt.savefig(os.path.join(os.path.dirname(paths['phaser']['imputed']), 'imputation_quality.png'))
plt.savefig(os.path.join(os.path.dirname(paths['beagle']['imputed']), 'imputation_quality.png'))
plt.show()
