import os
from typing import  *
import numpy as np
from cyvcf2 import VCF
from itertools import starmap, repeat
import shutil

from src.VCFPooling.poolSNPs import parameters as prm
from src.VCFPooling.poolSNPs import beagle_tools as bgltls

from persotools.files import delete_file, mkdir
from persotools.struct import NamedDict

'''
Run bcftools manipulations for preprocessing the files which are used for Beagle imputation.
Run Beagle.

v1.1: run imputation for single sample, one by one. (Compare with imputation on all together)
v1.1.1: implement multiprocessing

Steps:
* (transfer codes: different script)
* (transfer raw and pooled files to the server (adaptive GL schemes only): different script)
* split into reference panel (REF) and study population (IMP)
* split IMP into individual files
* phase and imput every individual file
* merge imputed files
* (transfer to local computer: different script)
'''

### PARAMETERS
chk_sz = prm.CHK_SZ
path_gt_files = prm.PATH_GT_FILES
path_gl_files = prm.PATH_GL_FILES

if prm.GTGL == 'GT':
    cd = path_gt_files
if prm.GTGL == 'GL':
    if prm.unknown_gl != 'adaptive':
        cd = os.path.join(path_gl_files, 'gl_' + '_'.join(np.around(prm.unknown_gl, decimals=2).astype(str)))
    else:
        cd = os.path.join(path_gl_files, 'gl_adaptive')
    mkdir(cd)
os.chdir(cd)

raw = NamedDict('raw', list(prm.RAW.keys()), list(prm.RAW.values()))
pooled = NamedDict('pooled', list(prm.POOLED.keys()), list(prm.POOLED.values()))

print('Start processing files for running BEAGLE'.ljust(80, '.'))
"""
### BGZIP ALL
bgltls.bgzip_working_files(pooled, path_gt_files, path_gl_files, cd)

### REF/IMP SAMPLING
bgltls.create_ref_imp_lists(cd, sizeref=prm.NB_REF, sizeimp=prm.NB_IMP)

for (folder, dic) in tuple([(cd, pooled),
                            (path_gt_files, raw)]):
    bgltls.partition_imp((folder, dic), total_ref=False)

bgltls.partition_ref(raw, path_gt_files)

for (folder, f) in tuple([(path_gt_files, raw['imp']),
                          (path_gt_files, raw['ref']),
                          (cd, pooled['imp'])]):
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
"""
"""
Important note!

In python 3.x, map creates an iterable. Iterables are evaluated lazily - that is,
only the elements you iterate over (such as looping over them) are actually evaluated.
If don't iterate over any values, none of the values get evaluated.
In the any case, however, any iterates over the entire result map, so the whole thing gets evaluated.

Alternative: queue module?
"""
### REFINE SAMPLES OR MARKERS TO IMPUTE
idv_to_keep: List[str] = VCF(os.path.join(path_gt_files, raw['imp'])).samples[0:3]
# .samples[0:1] gets only the first sample into a list i.e. ['NA18960']
# .samples[0] gets only the first sample as a string i.e. 'NA18960'
print(idv_to_keep)
args = list(zip(repeat(pooled, len(idv_to_keep)),
                repeat(cd, len(idv_to_keep)),
                idv_to_keep))
path_out = list(starmap(bgltls.keep_single_sample, args))
"""
### BEAGLE ROUND#1: (GL to GT and) PHASING
# TODO: bgltls.beagle_phasing(raw, path_gt_files, cd)
args = list(zip(repeat(pooled, len(path_out)),
                repeat(path_gt_files, len(path_out)),
                path_out))
_ = list(starmap(bgltls.beagle_phasing, args))

# for arg in args:  # REF-file access competition/restricted? use a queuing system with a simple for loop?
#    bgltls.beagle_phasing(*arg)

### CONFORM-GT
# conform-gt algo should be used if REF and IMP datasets contain different sets of markers
# chrom-POS: min and max from the original file
args = list(zip(repeat(pooled, len(path_out)),
                repeat(raw, len(path_out)),
                path_out))
_ = all(starmap(bgltls.conform_gt, args))

### BEAGLE (ROUND#2): IMPUTING
args = list(zip(repeat(pooled, len(path_out)),
                repeat(raw, len(path_out)),
                path_out))
_ = all(starmap(bgltls.beagle_imputing, args))
# for arg in args:
#     bgltls.beagle_imputing(*arg)  # _ = all(starmap(

# error: /bin/sh: 1: bgzip: not found
# bgzip is in the tabix package (not to confound with bzip2)
# sudo apt-get install tabix
# on Rackham: module load bioinfo-tools && module load tabix/0.2.6

### FIX DS AND GP FORMAT FIELDS, SORT OUT GT
args = list(zip(repeat(pooled, len(path_out)),
                path_out))
_ = all(starmap(bgltls.reformat_fields, args))

_ = all(map(bgltls.clean_imputed_directory, path_out))
_ = bgltls.clean_imputed_directory(cd)
"""
### MERGE UNIT FILES TOGETHER AND CLEAN DIRECTORIES
mkdir(os.path.join(cd, 'single_samples_merged'))
files2merge = []
for n in idv_to_keep:
    filein = os.path.join(cd,
                            'keeponly_{}'.format(n),
                            'IMP.chr20.{}.imputed.gt.chunk10000.vcf.gz'.format(pooled.name)
                            )
    fileout = os.path.join(cd, 'single_samples_merged', '{}.imputed.gt.unitfile.vcf.gz'.format(n))
    files2merge.append(fileout)
    _ = bgltls.move_file_to(filein, fileout, os.path.join(cd, 'single_samples_merged'))
    shutil.rmtree(os.path.join(cd, 'keeponly_{}'.format(n)))
_ = bgltls.merge_files('*.imputed.gt.unitfile.vcf.gz',
                       'IMP.chr20.{}.imputed.gt.chunk10000.vcf.gz'.format(pooled.name),
                       cd)
delete_file(os.path.join(cd, 'single_samples_merged', '*.imputed.gt.unitfile.vcf.gz'))
