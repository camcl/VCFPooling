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

v1.3: run imputation for all samples and all markers together with command-line launching. 
To be run on a remote server.

Calling from command-line: $ python3 beagle_impute_v1.3.py chr20_20200307 4

Note regarding threading with Beagle (from 4.1 documentation):
nthreads=[positive integer] specifies the number of threads of execution.  If no nthreads parameter  is  specified,
the nthreads  parameter  will  be  set  equal  to  the  number  of  CPU cores on the host machine. 

Steps:
* (transfer codes: other script to run locally)
* (transfer raw and pooled files to the server (adaptive GL schemes only): different script)
* split into reference panel (REF) and study population (IMP)
* split IMP into individual files
* phase and imput every individual file
* merge imputed files
* (transfer merged imputed files to local computer: other script to run locally)
'''
prm.info()

### COMMAND-LINE PARSING AND PARAMETERS
parser = argparse.ArgumentParser(description='Run imputation with Beagle for single sample,'
                                             'one by one and parallelized.')
parser.add_argument('cores', metavar='nc', type=int, nargs='?', help='Number of cores to use', default=None)
parser.add_argument('fpath', metavar='pth', type=str, help='Location of file to be imputed', default=None)

argsin = parser.parse_args()

nb_cores = os.cpu_count() if argsin.cores is None else argsin.cores
fpath = argsin.fpath

print('\n'.ljust(80, '*'))
print('Number of cpu to be used = {}'.format(nb_cores))
print('Process files in {}'.format(fpath))
print('\n'.rjust(80, '*'))

path_gt_files = prm.PATH_GT_FILES
path_gl_files = prm.PATH_GL_FILES

if prm.GTGL == 'GT':
    cd = path_gt_files
if prm.GTGL == 'GL':
    if prm.unknown_gl != 'adaptive':
        cd = os.path.join(path_gl_files, 'gl_' + '_'.join(np.around(prm.unknown_gl, decimals=2).astype(str)))
    else:
        cd = os.path.join(path_gl_files, 'gl_adaptive', fpath)
    mkdir(cd)
os.chdir(cd)

raw = NamedDict('raw', list(prm.RAW.keys()), list(prm.RAW.values()))
pooled = NamedDict('pooled', list(prm.POOLED.keys()), list(prm.POOLED.values()))

print('\nStart processing files for running BEAGLE'.ljust(80, '.'))

### BGZIP ALL
# bgltls.bgzip_working_files(pooled, path_gt_files, path_gl_files, cd)

### REF/IMP SAMPLING
bgltls.create_ref_imp_lists(cd, sizeref=prm.NB_REF, sizeimp=prm.NB_IMP)

for (folder, dic) in tuple([(cd, pooled),
                            (path_gt_files, raw)]):
    bgltls.partition_imp((folder, dic), total_ref=False)

#bgltls.partition_ref(raw, path_gt_files)
bgltls.partition_ref(raw, cd)

for (folder, f) in tuple([(path_gt_files, raw['imp']),
                          (path_gt_files, raw['ref']),
                          (cd, pooled['imp'])]):
    pybcf.sort(f, folder)
    pybcf.index(f, folder)

if prm.GTGL == 'GT':
    bgltls.partition_ref(pooled, path_gt_files)

"""
Important note!

In python 3.x, map creates an iterable. Iterables are evaluated lazily - that is,
only the elements you iterate over (such as looping over them) are actually evaluated.
--> Use all(iterable) or list(iterable) to evaluate 
"""
### REFINE SAMPLES OR MARKERS TO IMPUTE
idv_to_keep: List[str] = VCF(os.path.join(path_gt_files, raw['imp'])).samples
print('List of samples to be imputed --> {}'.format(idv_to_keep))
path_out = cd

### BEAGLE ROUND#1: (GL to GT and) PHASING
_ = bgltls.beagle_phasing(raw, path_gt_files, path_out)
_ = bgltls.beagle_phasing(pooled, path_gt_files, path_out)

### CONFORM-GT
# conform-gt algo should be used if REF and IMP datasets contain different sets of markers
# chrom-POS: min and max from the original file
_ = bgltls.conform_gt(pooled, raw, path_out)

### BEAGLE (ROUND#2): IMPUTING
_ = bgltls.beagle_imputing(pooled, raw, path_out)

# error: /bin/sh: 1: bgzip: not found
# bgzip is in the tabix package (not to confound with bzip2)
# sudo apt-get install tabix
# on Rackham: module load bioinfo-tools && module load tabix/0.2.6

### FIX DS AND GP FORMAT FIELDS, SORT OUT GT
_ = bgltls.reformat_fields(pooled, path_out)
_ = bgltls.clean_imputed_directory(path_out)
_ = bgltls.clean_imputed_directory(cd)

### MERGE UNIT FILES TOGETHER AND CLEAN DIRECTORIES
pass
