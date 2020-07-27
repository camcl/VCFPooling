import sys, os
sys.path.insert(0, '/home/camille/1000Genomes/src/')
import pandas as pd
import numpy as np

from VCFPooling.poolSNPs.metrics import quality
import subprocess
import argparse
import matplotlib.pyplot as plt


"""
Compute results with customized metrics from true vs. imputed data sets
ex.
$ python3 -u imputation_quality.py<path to directory> <VCF file with true genotypes> <VCF file with imputed genotypes> <path to script converting GT to GL>
"""

parser = argparse.ArgumentParser(description='Compute and plot'
                                             'customized imputation accuracy metrics')
parser.add_argument('directory', metavar='dir', type=str, help='Path to directory with files', default=None)
parser.add_argument('true', metavar='tru', type=str, help='File with true genotypes', default=None)
parser.add_argument('imputed', metavar='imp', type=str, help='Imputed file with genotypes (GT:DS:GP)', default=None)
parser.add_argument('gconverter', metavar='gcv', type=str, help='Imputed file with genotypes (GT:DS:GP)', default='~/PoolImpHuman/bin/bash-src/gt_to_gl.sh')

argsin = parser.parse_args()
dirin = argsin.directory
ftrue = argsin.true
fimp = argsin.imputed
gconv = argsin.gconverter

paths = {'beaglegt': {
    'true': os.path.join(dirin, ftrue),
    'imputed': os.path.join(dirin, fimp)},
    'beaglegl': {
        'true': os.path.join(dirin, ftrue.replace('.gt.', '.gl.')),
        'imputed': os.path.join(dirin, fimp)},
}

convertgtgl = True
if convertgtgl:
    cmd = 'bash {} {} {}'.format(gconv, paths['beaglegt']['true'], paths['beaglegl']['true'])
    subprocess.run(cmd, shell=True,)

qbeaglegt = quality.QualityGT(*paths['beaglegt'].values(), 0, idx='id')
qbeaglegl = quality.QualityGL(paths['beaglegl']['true'], paths['beaglegl']['imputed'], 0, idx='id')
messbeagle = qbeaglegl.cross_entropy

tabbeaglegl = pd.concat([qbeaglegt.concordance(),
                         qbeaglegt.trueobj.af_info,
                         qbeaglegt.pearsoncorrelation(),
                         qbeaglegt.precision,
                         qbeaglegt.accuracy,
                         qbeaglegt.recall,
                         qbeaglegt.f1_score,
                         qbeaglegl.cross_entropy], axis=1)
dosbeaglegl = qbeaglegt.alleledosage()

tabbeaglegl.head()

plt.rcParams["figure.figsize"] = [5*4, 4*2]
fig, axes = plt.subplots(2, 4)

tabbeaglegl.plot.scatter('af_info', 'precision_score', ax=axes[0, 0], s=0.7, label='beaglegl')
axes[0, 0].set_ylim(0.0, 1.0)
tabbeaglegl.plot.scatter('af_info', 'accuracy_score', ax=axes[0, 1], s=0.7, label='beaglegl')
axes[0, 1].set_ylim(0.0, 1.0)
tabbeaglegl.plot.scatter('af_info', 'concordance', ax=axes[0, 2], s=0.7, label='beaglegl')
axes[0, 2].set_ylim(0.0, 1.0)
tabbeaglegl.plot.scatter('af_info', 'f1_score', ax=axes[0, 3], s=0.7, label='beaglegl')
axes[0, 3].set_ylim(0.0, 1.0)
tabbeaglegl.plot.scatter('af_info', 'r_squared', ax=axes[1, 0], s=0.7, label='beaglegl')
axes[1, 0].set_ylim(-0.2, 1.0)
tabbeaglegl.plot.scatter('af_info', 'cross_entropy', ax=axes[1, 1], s=0.7, label='beaglegl')
axes[1, 1].set_ylim(-0.5, 5.0)
axes[1, 2].scatter(dosbeaglegl[0], dosbeaglegl[1], s=0.7, label='beaglegl')
axes[1, 2].set_xlabel('true allele dosage')
axes[1, 2].set_ylabel('imputed allele dosage')
axes[1, 2].set_ylim(0.0, 2.0)

for ax in axes.flatten()[:-2]:
    # cast color to white 'w' if dark background
    ax.set_xlabel('true alternate allele frequency', color='w')
    ax.set_ylabel(ax.get_ylabel(), color='w')
plt.savefig(os.path.join(os.path.dirname(paths['beaglegt']['imputed']), 'imputation_quality_gtgl.png'))
plt.show()
