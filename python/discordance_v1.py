import numpy as np
import itertools

from scripts.VCFPooling.python.archived.alleles import alleles_tools as alltls
from scripts.VCFPooling.poolSNPs import dataframe as vcfdf
from scripts.VCFPooling.poolSNPs import parameters as prm
from persotools.files import *

"""
Compute discordance rates for different data sets imputed.

Data sets:
0) GL settings from 201908
To be compared:
1) pooled data imputed from adaptive GL with Beagle, all samples from study pop imputed together
vs.
2) pooled data imputed from adaptive GL with Beagle, samples from study pop imputed one-by-one
"""

print('Configure paths'.ljust(80, '.'))
path_gt_beagle = os.path.join(prm.WD, 'gt', 'stratified', 'all_snps_all_samples')
path_gl_beagle1 = os.path.join(prm.WD, 'gl', 'gl_adaptive', 'chunk10000_20190822')
path_gl_beagle2 = os.path.join(prm.WD, 'gl', 'gl_adaptive', 'single_samples_merged')

paths = [path_gt_beagle,
         path_gl_beagle1,
         path_gl_beagle2]

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
pdvcf = vcfdf.PandasVCF(os.path.join(prm.PATH_GT_FILES, prm.CHKFILE), indextype='id')
aafs = pdvcf.concatcols([pdvcf.af_info, pdvcf.aaf])

# Create line-index and column-index (multi-indexing)
print('Create multi-indices for the study population (IMP data sets)'.ljust(80, '.'))
aaf_idx, pop_idx = alltls.make_index(RAW)

# Load data sets for each working directory
vcfraw = vcfdf.PandasVCF(RAW, indextype='id')
raw0, raw1 = vcfraw.vcf2dframe()
samples = vcfraw.samples
variants = vcfraw.variants
raw = np.stack([raw0.values, raw1.values], axis=-1)

for cd in vcfpathdic.keys():
    os.chdir(cd)
    print('\r\nBuild datasets in {}'.format(cd).ljust(80, '.'))
    # Raw data set as ground truth for genotypes in the study population
    for name, dset in zip(setnames, vcfpathdic[cd]):
        if dset is not None:
            print('Error and discordance in {}'.format(dset).ljust(80, '.'))
            vcfset = vcfdf.PandasVCF(dset)
            dset0, dset1 = vcfset.vcf2dframe()

            # NumPy juxtapos for error computation
            dset = np.stack([dset0.values, dset1.values], axis=-1)

            set_discordance = alltls.compute_discordance([name],
                                                         [dset],
                                                         raw,
                                                         aaf_idx[:raw.shape[0]],
                                                         pop_idx,
                                                         sorting)
            print('Errors from alltls.compute_discordance --> ', set_discordance[name].mean().mean())


