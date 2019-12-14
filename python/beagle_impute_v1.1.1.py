import sys, os
from typing import *
import numpy as np
from cyvcf2 import VCF
from itertools import starmap, repeat
import shutil
import multiprocessing as mp
import argparse

home_dir = os.path.expanduser("~")
proj_dir = os.path.join(home_dir, '1000Genomes')
sys.path.insert(0, proj_dir)

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
v1.1.1: implement multiprocessing and command-line launching. To be run from a remote server

Calling from command-line: $ python3 beagle_impute_v1.1.1.py 4 8

Steps:
* (transfer codes: other script to run locally)
* (transfer raw and pooled files to the server (adaptive GL schemes only): different script)
* split into reference panel (REF) and study population (IMP)
* split IMP into individual files
* phase and imput every individual file
* merge imputed files
* (transfer merged imputed files to local computer: other script to run locally)
'''

### COMMAND-LINE PARSING AND PARAMETERS
parser = argparse.ArgumentParser(description='Run imputation with Beagle for single sample,'
                                             'one by one and parallelized.')
parser.add_argument('cores', metavar='nc', type=int, nargs='?', help='Number of cores to use', default=None)
parser.add_argument('idv_nb', metavar='ni', type=int, nargs='?', help='Number of samples to impute', default=None)
argsin = parser.parse_args()

nb_cores = os.cpu_count() if argsin.cores is None else argsin.cores
nb_idv = prm.NB_IMP if argsin.idv_nb is None else argsin.idv_nb

print('\n'.ljust(80, '*'))
print('Number of cpu to be used = {}'.format(nb_cores))
print('Number of samples to be imputed = {}'.format(nb_idv))
print('\n'.rjust(80, '*'))

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

print('\nStart processing files for running BEAGLE'.ljust(80, '.'))
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
--> Use all(iterable) or list(iterable) to evaluate 
"""
### REFINE SAMPLES OR MARKERS TO IMPUTE
idv_to_keep: List[str] = VCF(os.path.join(path_gt_files, raw['imp'])).samples[0:nb_idv]
print(idv_to_keep)
args = list(zip(repeat(pooled, len(idv_to_keep)),
                repeat(cd, len(idv_to_keep)),
                idv_to_keep))
with mp.Pool(processes=nb_cores) as mpool:
    path_out = list(mpool.starmap(bgltls.keep_single_sample, args))
### BEAGLE ROUND#1: (GL to GT and) PHASING
# TODO: bgltls.beagle_phasing(raw, path_gt_files, cd)
args = list(zip(repeat(pooled, len(path_out)),
                repeat(path_gt_files, len(path_out)),
                path_out))
with mp.Pool(processes=nb_cores) as mpool:
    _ = list(mpool.starmap(bgltls.beagle_phasing, args))

### CONFORM-GT
# conform-gt algo should be used if REF and IMP datasets contain different sets of markers
# chrom-POS: min and max from the original file
args = list(zip(repeat(pooled, len(path_out)),
                repeat(raw, len(path_out)),
                path_out))
with mp.Pool(processes=nb_cores) as mpool:
    _ = all(mpool.starmap(bgltls.conform_gt, args))

### BEAGLE (ROUND#2): IMPUTING
args = list(zip(repeat(pooled, len(path_out)),
                repeat(raw, len(path_out)),
                path_out))
with mp.Pool(processes=nb_cores) as mpool:
    _ = all(mpool.starmap(bgltls.beagle_imputing, args))

# error: /bin/sh: 1: bgzip: not found
# bgzip is in the tabix package (not to confound with bzip2)
# sudo apt-get install tabix
# on Rackham: module load bioinfo-tools && module load tabix/0.2.6

### FIX DS AND GP FORMAT FIELDS, SORT OUT GT
args = list(zip(repeat(pooled, len(path_out)),
                path_out))
with mp.Pool(processes=nb_cores) as mpool:
    _ = all(mpool.starmap(bgltls.reformat_fields, args))
    _ = all(mpool.map(bgltls.clean_imputed_directory, path_out))
_ = bgltls.clean_imputed_directory(cd)

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
