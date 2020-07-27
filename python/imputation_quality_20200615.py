import matplotlib.pyplot as plt
import pandas as pd
import subprocess

from src.VCFPooling.poolSNPs.metrics import quality

from persotools.files import *

"""
Compare imputation results from Beagle on pooled file obtained from the Pysam-based implementation.
Plot imputation performance in:
- precision
- accuracy
- recall,
- f1_score
- allele dosage
- Pearson's correlation
"""

convert_files = True

# Example with results from Beagle on 10000 markers
paths = {'beagle1': {
    'true': '/home/camille/PoolImpHuman/data/20200615/IMP.chr20.snps.gt.chunk10000.vcf.gz',
    'imputed': '/home/camille/PoolImpHuman/data/20200615/IMP.chr20.pooled.imputed.gt.chunk10000.vcf.gz'},
         'beaglegl': {
     'true': '/home/camille/PoolImpHuman/data/20200615/IMP.chr20.snps.gl.chunk10000.vcf.gz',
     'imputed': '/home/camille/PoolImpHuman/data/20200615/IMP.chr20.pooled.beagle2.gl.chunk10000.vcf.gz'},
}
# Convert GT files to GL (Phaser)
if convert_files:
    cmd = 'bash ~/PoolImpHuman/bin/bash-src/gt_to_gl.sh {} {}'.format(paths['beagle1']['true'], paths['beaglegl']['true'])
    subprocess.run(cmd, shell=True,)
    # easier for cross entropies not to log gl

qbeagle1 = quality.QualityGT(*paths['beagle1'].values(), 0, idx='id')
tabbeagle1 = pd.concat([qbeagle1.concordance(),
                       qbeagle1.trueobj.af_info,
                       qbeagle1.pearsoncorrelation(),
                       qbeagle1.precision,
                       qbeagle1.accuracy,
                       qbeagle1.recall,
                       qbeagle1.f1_score], axis=1)
dosbeagle1 = qbeagle1.alleledosage()

qbeaglegl = quality.QualityGL(paths['beaglegl']['true'], paths['beaglegl']['imputed'], 0, idx='id')
messbeagle = qbeaglegl.cross_entropy

qbeagle2 = quality.QualityGT(*paths['beagle1'].values(), 0, idx='id')
tabbeagle2 = pd.concat([qbeagle2.concordance(),
                       qbeagle2.trueobj.af_info,
                       qbeagle2.pearsoncorrelation(),
                       qbeagle2.precision,
                       qbeagle2.accuracy,
                       qbeagle2.recall,
                       qbeagle2.f1_score,
                        qbeaglegl.cross_entropy], axis=1)
dosbeagle2 = qbeagle2.alleledosage()

plt.rcParams["figure.figsize"] = [5*2, 4*7]
fig, axes = plt.subplots(7, 2)

tabbeagle2.plot.scatter('af_info', 'precision_score', ax=axes[0, 0], s=0.7, label='beagle2')
axes[0, 0].set_ylim(0.0, 1.0)
tabbeagle2.plot.scatter('af_info', 'accuracy_score', ax=axes[1, 0], s=0.7, label='beagle2')
axes[1, 0].set_ylim(0.0, 1.0)
tabbeagle2.plot.scatter('af_info', 'concordance', ax=axes[2, 0], s=0.7, label='beagle2')
axes[2, 0].set_ylim(0.0, 1.6)
tabbeagle2.plot.scatter('af_info', 'f1_score', ax=axes[3, 0], s=0.7, label='beagle2')
axes[3, 0].set_ylim(0.0, 1.0)
tabbeagle2.plot.scatter('af_info', 'r_squared', ax=axes[4, 0], s=0.7, label='beagle2')
axes[4, 0].set_ylim(-0.4, 1.0)
tabbeagle2.plot.scatter('af_info', 'cross_entropy', ax=axes[5, 0], s=0.7, label='beagle2')
axes[5, 0].set_ylim(0.0, 1.2)
axes[6, 0].scatter(dosbeagle2[0], dosbeagle2[1], s=0.7, label='beagle2')
axes[6, 0].set_xlabel('true allele dosage')
axes[6, 0].set_ylabel('imputed allele dosage')
axes[6, 0].set_ylim(0.0, 2.0)

tabbeagle1.plot.scatter('af_info', 'precision_score', ax=axes[0, 1], s=0.7, label='beagle1')
axes[0, 1].set_ylim(0.0, 1.0)
tabbeagle1.plot.scatter('af_info', 'accuracy_score', ax=axes[1, 1], s=0.7, label='beagle1')
axes[1, 1].set_ylim(0.0, 1.0)
tabbeagle1.plot.scatter('af_info', 'concordance', ax=axes[2, 1], s=0.7, label='beagle1')
axes[2, 1].set_ylim(0.0, 1.6)
tabbeagle1.plot.scatter('af_info', 'f1_score', ax=axes[3, 1], s=0.7, label='beagle1')
axes[3, 1].set_ylim(0.0, 1.0)
tabbeagle1.plot.scatter('af_info', 'r_squared', ax=axes[4, 1], s=0.7, label='beagle1')
axes[4, 1].set_ylim(-0.4, 1.0)
axes[5, 1].set_ylim(0.0, 1.2)
axes[6, 1].scatter(dosbeagle1[0], dosbeagle1[1], s=0.7, label='beagle1')
axes[6, 1].set_xlabel('true allele dosage')
axes[6, 1].set_ylabel('imputed allele dosage')
axes[6, 1].set_ylim(0.0, 2.0)
for ax in axes.flatten()[:-2]:
    ax.set_xlabel('true alternate allele frequency')
# plt.savefig(os.path.join(os.path.dirname(paths['beagle1']['imputed']), 'imputation_quality.png'))
plt.savefig(os.path.join(os.path.dirname(paths['beagle1']['imputed']), 'imputation_quality_gl.png'))
plt.show()
