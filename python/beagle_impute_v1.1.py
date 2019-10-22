import os
from typing import  *
import numpy as np
from cyvcf2 import VCF
from itertools import starmap, repeat

from scripts.VCFPooling.poolSNPs import parameters as prm
from scripts.VCFPooling.poolSNPs import beagle_tools as bgltls
from scripts.VCFPooling.poolSNPs import pybcf
from scripts.VCFPooling.poolSNPs.alleles import alleles_tools as alltls

from persotools.files import delete_file, mkdir
from persotools.struct import NamedDict

'''
Run bcftools manipulations for preprocessing the files which are used for Beagle imputation.
Run Beagle.

v1.1: run imputation for single sample, one by one. (Compare with imputation on all together)
v1.1.1: implement multiprocessing
'''

### PARAMETERS
chk_sz = prm.CHK_SZ
path_gt_files = prm.PATH_GT_FILES
path_gl_files = prm.PATH_GL_FILES

if prm.GTGL == 'GT':
    cd = path_gt_files
if prm.GTGL == 'GL':
    if prm.unknown_gl != 'adaptative':
        cd = os.path.join(path_gl_files, 'gl_' + '_'.join(np.around(prm.unknown_gl, decimals=2).astype(str)))
    else:
        cd = os.path.join(path_gl_files, 'gl_adaptative')
    mkdir(cd)
os.chdir(cd)

raw = NamedDict('raw', list(prm.RAW.keys()), list(prm.RAW.values()))
pooled = NamedDict('pooled', list(prm.POOLED.keys()), list(prm.POOLED.values()))
missing = NamedDict('missing', list(prm.MISSING.keys()), list(prm.MISSING.values()))
subset = prm.SUBSET #if susbset get the first 1000 lines to recreate ech method file: pooled, missing etc
trc = prm.SUBCHUNK

print('Start processing files for running BEAGLE'.ljust(80, '.'))

if subset:
    print('Subset main set in {}'.format(os.getcwd()).ljust(80, '.'))
    for dic in [pooled, missing, raw]:
        for k, v in dic.items():
            dic[k] = v.replace('chunk' + str(chk_sz), 'chunk' + str(trc))

### BGZIP ALL
for dic in [pooled, missing]:
    bgltls.bgzip_working_files(dic, path_gt_files, path_gl_files, cd)

### REF/IMP SAMPLING
bgltls.create_ref_imp_lists(cd, sizeref=prm.NB_REF, sizeimp=prm.NB_IMP)

for (folder, dic) in tuple([(cd, pooled),
                            (cd, missing),
                            (path_gt_files, raw)]):
    bgltls.partition_imp((folder, dic), total_ref=False)

bgltls.partition_ref(raw, path_gt_files)

for (folder, f) in tuple([(path_gt_files, raw['imp']),
                          (path_gt_files, raw['ref']),
                          (cd, pooled['imp']),
                          (cd, missing['imp'])]):
    pybcf.sort(f, folder)
    pybcf.index(f, folder)

if prm.GTGL == 'GT':
    bgltls.partition_ref(pooled, path_gt_files)

### GL CONVERSION FOR IMP/REF
if prm.GTGL == 'G':
    print('GT to GL in {}'.format(os.getcwd()).ljust(80, '.'))
    alltls.file_likelihood_converter(os.path.join(path_gt_files,
                                                  raw['ref'].replace('.gl', '.gt')),
                                     raw['ref'][:-3])
    delete_file(raw['ref'])
    delete_file(raw['ref'] + '.csi')
    pybcf.bgzip(raw['ref'][:-3], raw['ref'], cd)
    pybcf.sort(raw['ref'], cd)
    pybcf.index(raw['ref'], cd)
    delete_file(raw['ref'][:-3])

### BEAGLE ROUND#1: (GL to GT and) PHASING
for dic in [raw, pooled, missing]:
    bgltls.beagle_haplo_to_geno()
    bgltls.beagle_phasing(dic, path_gt_files, cd)

"""
Important note!

In python 3.x, map creates an iterable. Iterables are evaluated lazily - that is,
only the elements you iterate over (such as looping over them) are actually evaluated.
If don't iterate over any values, none of the values get evaluated.
In the any case, however, any iterates over the entire result map, so the whole thing gets evaluated.

Alternative: queue module
"""
### REFINE SAMPLES OR MARKERS TO IMPUTE
snp_to_rm: List[str] = []
idv_to_keep: List[str] = VCF(os.path.join(path_gt_files, raw['imp'])).samples

for dic in [pooled, missing]:
    if len(snp_to_rm) > 0:
        args = list(zip(repeat(dic, len(snp_to_rm)),
                        repeat(cd, len(snp_to_rm)),
                        snp_to_rm))
        path_out = list(starmap(bgltls.switch_off_markers, args))
    elif len(idv_to_keep) > 0:
        args = list(zip(repeat(dic, len(idv_to_keep)),
                        repeat(cd, len(idv_to_keep)),
                        idv_to_keep))
        path_out = list(starmap(bgltls.keep_single_sample, args))
    else:
        path_out = bgltls.all_snps_all_samples(dic, cd)

### CONFORM-GT
# conform-gt algo should be used if REF and IMP datasets contain different sets of markers
# chrom-POS: min and max from the original file
for dic in [pooled, missing]:
    args = list(zip(repeat(dic, len(path_out)),
                    repeat(raw, len(path_out)),
                    path_out))
    _ = all(starmap(bgltls.conform_gt, args))

### BEAGLE (ROUND#2): IMPUTING
for dic in [pooled, missing]:
    args = list(zip(repeat(dic, len(path_out)),
                    repeat(raw, len(path_out)),
                    path_out))
    _ = all(starmap(bgltls.beagle_imputing, args))

### FIX DS AND GP FORMAT FIELDS, SORT OUT GT
for dic in [pooled, missing]:
    args = list(zip(repeat(dic, len(path_out)),
                    path_out))
    _ = all(starmap(bgltls.reformat_fields, args))

_ = all(map(bgltls.clean_imputed_directory, path_out))
_ = bgltls.clean_imputed_directory(cd)

### MERGE UNIT FILES TOGETHER
for dic in [pooled, missing]:
    files2merge = list([os.path.join('../keeponly_{}'.format(n),
                                     'IMP.chr20.{}.imputed.gt.chunk10000.vcf.gz'.format(dic.name)
                                     )
                        for n in idv_to_keep])
    _ = bgltls.merge_files(files2merge,
                           'IMP.chr20.{}.imputed.gt.chunk10000.vcf.gz'.format(dic.name),
                           cd)