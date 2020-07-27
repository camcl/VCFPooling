import matplotlib.pyplot as plt
import pandas as pd

from src.VCFPooling.poolSNPs.metrics import quality

from persotools.files import *

"""
Plot imputation performance in:
- precision
- accuracy
- recall,
- f1_score
- allele dosage
- Pearson's correlation
"""

# Example with results from Phaser and Beagle on 10000 markers
paths = {'phaser': {
    'true': '/home/camille/1000Genomes/data/gt/stratified/IMP.chr20.snps.gt.chunk10000.vcf.gz',
    'imputed': '/home/camille/1000Genomes/data/gl/gl_adaptive/phaser/IMP.chr20.pooled.imputed.gt.chunk10000.vcf.gz'},
         'beagle': {
     'true': '/home/camille/1000Genomes/data/gt/stratified/IMP.chr20.snps.gt.chunk10000.vcf.gz',
     'imputed': '/home/camille/1000Genomes/data/gl/gl_adaptive/chunk10000_20190725/IMP.chr20.pooled.imputed.gt.chunk10000.vcf.gz'}
}
qphaser = quality.QualityGT(*paths['phaser'].values(), 0, idx='chrom:pos')
tabphaser = pd.concat([qphaser.concordance(),
                       qphaser.trueobj.af_info,
                       qphaser.pearsoncorrelation(),
                       qphaser.precision,
                       qphaser.accuracy,
                       qphaser.recall,
                       qphaser.f1_score], axis=1)
dosphaser = qphaser.alleledosage()

qbeagle = quality.QualityGT(*paths['beagle'].values(), 0, idx='chrom:pos')
tabbeagle = pd.concat([qbeagle.concordance(),
                       qbeagle.trueobj.af_info,
                       qbeagle.pearsoncorrelation(),
                       qbeagle.precision,
                       qbeagle.accuracy,
                       qbeagle.recall,
                       qbeagle.f1_score], axis=1)
dosbeagle = qbeagle.alleledosage()

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
# tabbeagle.plot.scatter('af_info', 'cross_entropy', ax=axes[5, 0], s=0.7, label='beagle')
# axes[5, 0].set_ylim(0.0, 1.1)
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
# tabphaser.plot.scatter('af_info', 'cross_entropy', ax=axes[5, 1], s=0.7, label='phaser')
# axes[5, 1].set_ylim(0.0, 1.1)
axes[6, 1].scatter(dosphaser[0], dosphaser[1], s=0.7, label='phaser')
axes[6, 1].set_xlabel('true allele dosage')
axes[6, 1].set_ylabel('imputed allele dosage')
axes[6, 1].set_ylim(0.0, 2.0)
for ax in axes.flatten()[:-2]:
    ax.set_xlabel('true alternate allele frequency')
plt.savefig(os.path.join(os.path.dirname(paths['phaser']['imputed']), 'imputation_quality.png'))
plt.savefig(os.path.join(os.path.dirname(paths['beagle']['imputed']), 'imputation_quality.png'))
plt.show()
