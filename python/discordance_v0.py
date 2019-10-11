import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from cyvcf2 import VCF
import itertools

from scripts.poolSNPs.alleles import alleles_tools as alltls
from scripts.poolSNPs import parameters as prm
from scripts.poolSNPs import utils
from persotools.files import *

"""
Compute discordance rates for different data sets imputed.
Resulkts to be presented first at the Essence 2019 conference in Lund.

Data sets:
1) randomly missing data imputed from GT with Beagle
2) pooled data imputed from GT with Beagle
3) pooled data imputed from adaptive GL with Beagle
4) pooled data imputed from adaptive GL with Phaser
"""

print('Configure paths'.ljust(80, '.'))
path_gt_beagle = os.path.join(prm.WD, 'gt', 'stratified', 'all_snps_all_samples')
path_gl_beagle = os.path.join(prm.WD, 'gl', 'gl_adaptive', 'chunk10000_20190725')
# path_gl_beagle = os.path.join(prm.WD, 'gl', 'gl_adaptive', 'all_snps_all_samples')
path_gl_phaser = os.path.join(prm.WD, 'gl', 'gl_adaptive', 'phaser')

paths = [path_gt_beagle,
         path_gl_beagle,
         path_gl_phaser]

ALL = os.path.join(prm.WD, 'gt', 'stratified', prm.CHKFILE)
RAW = os.path.join(prm.WD, 'gt', 'stratified', prm.RAW['b1i'] + '.vcf.gz')
MISS = [os.path.join(path_gt_beagle, prm.MISSING['gtonly'] + '.vcf.gz')]
POOL = []
for cd in paths:
    POOL.append(os.path.join(cd, prm.POOLED['gtonly'] + '.vcf.gz'))

# Group data sets to process for each working directory
vcfpathdic = dict(zip(paths, itertools.zip_longest(POOL, MISS, fillvalue=None)))
setnames = ['pooled', 'missing']
print('Directories and data sets to be processed --> {}'.format(vcfpathdic))

print('Configure parameters'.ljust(80, '.'))
sorting = True
popsize = prm.NB_IMP
chk_sz = prm.CHK_SZ

# Load AAFs
print('Load aaf from {}'.format(ALL).ljust(80, '.'))
aafs = alltls.get_aaf(prm.WD + '/gt/stratified/ALL.chr20.snps.gt.chunk{}.vcf.gz'.format(str(chk_sz)))

# Create line-index and column-index (multi-indexing)
print('Create multi-indices for the study population (IMP data sets)'.ljust(80, '.'))
aaf_idx, pop_idx = alltls.make_index(RAW)

# Load AAFs
print('Load alleles frequencies'.ljust(80, '.'))
aafs = alltls.get_aaf(prm.WD + '/gt/stratified/ALL.chr20.snps.gt.chunk{}.vcf.gz'.format(str(chk_sz)))

# Load data sets for each working directory
raw0, raw1 = alltls.vcf2dframe(VCF(RAW), popsize)
raw0, raw1 = utils.sort_datasets([raw0, raw1], pop_idx, aafs)
samples = VCF(RAW).samples
variants = raw0.index.tolist()
raw = np.stack([raw0.values, raw1.values], axis=-1)

for cd in vcfpathdic.keys():
    os.chdir(cd)
    print('\r\nBuild datasets in {}'.format(cd).ljust(80, '.'))
    # Raw data set as ground truth for genotypes in the study population
    for name, dset in zip(setnames, vcfpathdic[cd]):
        if dset is not None:
            print('Error and discordance in {}'.format(dset).ljust(80, '.'))
            dset0, dset1 = alltls.vcf2dframe(VCF(dset), popsize)
            dset0, dset1 = utils.sort_datasets([dset0, dset1], pop_idx, aafs)

            # NumPy juxtapos for error computation
            dset = np.stack([dset0.values, dset1.values], axis=-1)

            set_discordance = alltls.compute_discordance([name],
                                                         [dset],
                                                         raw,
                                                         aaf_idx[:raw.shape[0]],
                                                         pop_idx,
                                                         sorting)
            print('Errors from alltls.compute_discordance --> ', set_discordance[name].mean().mean())

