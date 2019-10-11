import os
from typing import *
import numpy as np
from itertools import repeat
import multiprocessing as mp

from scripts.poolSNPs import parameters as prm
from scripts.poolSNPs import beagle_tools as bgltls
from scripts.poolSNPs import pybcf
from scripts.poolSNPs.alleles import alleles_tools as alltls

from persotools.files import delete_file, mkdir
from persotools.struct import NamedDict

'''
Run bcftools manipulations for preprocessing the files which are used for Beagle imputation.
Run Beagle.

v1.2: run imputation for 1 marker deleted. (Compare with imputation on all together)
v1.2.1: implements multiprocessing. Preparatory refactoring for working on Uppmax servers.
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
subset = prm.SUBSET  # if susbset get the first 1000 lines to recreate ech method file: pooled, missing etc
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
nb_cores = os.cpu_count() // 2
print('System\'s CPU number = {}'.format(nb_cores))
candidates = alltls.extract_variant_onfreq(os.path.join(path_gt_files,
                                                        raw['ref'].replace('.gl', '.gt')),
                                           [0.0, 1.0])
snp_to_rm: List[str] = candidates['id'].str.strip(' ').str.slice(3).values[:4]
idv_to_keep: List[str] = []
print('Number of SNPs to analyze separately: ', len(snp_to_rm))  # 225 with [0.49, 0.51], 121 with [0.495, 0.505]
print(snp_to_rm)

for dic in [pooled, missing]:
    with mp.Pool(processes=nb_cores) as workers:
        if len(snp_to_rm) > 0:
            args = list(zip(repeat(dic, len(snp_to_rm)),
                            repeat(cd, len(snp_to_rm)),
                            snp_to_rm))
            path_out = list(workers.starmap(bgltls.switch_off_markers, args))
        elif len(idv_to_keep) > 0:
            args = list(zip(repeat(dic, len(idv_to_keep)),
                            repeat(cd, len(idv_to_keep)),
                            idv_to_keep))
            path_out = list(workers.starmap(bgltls.keep_single_sample, args))
        else:
            path_out = bgltls.all_snps_all_samples(dic, cd)

### CONFORM-GT
# conform-gt algo should be used if REF and IMP datasets contain different sets of markers
# chrom-POS: min and max from the original file
for dic in [pooled, missing]:
    with mp.Pool(processes=nb_cores) as workers:
        args = list(zip(repeat(dic, len(path_out)),
                    repeat(raw, len(path_out)),
                    path_out))
        _ = all(workers.starmap(bgltls.conform_gt, args))

### BEAGLE (ROUND#2): IMPUTING
for dic in [pooled, missing]:
    with mp.Pool(processes=nb_cores) as workers:
        args = list(zip(repeat(dic, len(path_out)),
                    repeat(raw, len(path_out)),
                    path_out))
        _ = all(workers.starmap(bgltls.beagle_imputing, args))

### FIX DS AND GP FORMAT FIELDS, SORT OUT GT
for dic in [pooled, missing]:
    with mp.Pool(processes=nb_cores) as workers:
        args = list(zip(repeat(dic, len(path_out)),
                    path_out))
        _ = all(workers.starmap(bgltls.reformat_fields, args))
        _ = all(workers.map(bgltls.clean_imputed_directory, path_out))
bgltls.clean_imputed_directory(cd)

### MERGE UNIT FILES TOGETHER
for dic in [pooled, missing]:
    files2merge = list([os.path.join(cd,
                                     'rm_20:{}'.format(n),
                                     'marker_20:{}.vcf.gz'.format(n))
                        for n in snp_to_rm])
    bgltls.concat_files(files2merge,
                        '{}.vcf.gz'.format(dic['gtonly']),
                        dic,
                        cd)
