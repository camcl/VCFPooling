#!/usr/bin/env python

#TODO: Rename files with core + extension (.vcf.gz, .vcf, .csv, ...)
#TODO: files_root as parameter and create the whole tree structure (independant from ~PycharmProject/1000Genomes
"""
Parameters for running pooling + decoding on a dataset:
- path to file to process
- make new directory if needed: today's date, split data and pictures?
- files names/extension for outputs
- iteration number missing/pooled/missing.pooled
- chunks size if needed
"""

import os
import numpy as np
from cyvcf2 import VCF
from typing import *
from pathlib import Path


### GENERAL TOOLS


def pgcd(a, b):
    """Computes the biggest common divider"""
    while b != 0:
        r = a % b
        a, b = b, r
    return a


def ppcm(a, b):
    """Computes the smallest common multiple"""
    if (a == 0) or (b == 0):
        return 0
    else:
        return (a*b)//pgcd(a, b)


ext_vcf = '.vcf'
ext_vcfgz = '.vcf.gz'
ext_csv = '.csv'

GLtype = NewType('GLtype', Tuple[float, float, float])
GTtype = NewType('GTtype', Tuple[int, int, bool])

### tree struct
HOME = str(Path.home())
# TODO: replace all occurrences of /home/camille by prm.HOME
ROOT = os.path.join(HOME, '1000Genomes')

### all
DATA_PATH = os.path.join(ROOT, 'data')
PLOTS_PATH = os.path.join(ROOT, 'plots')
SCRIPTS_PATH = os.path.join(ROOT, 'scripts')


### pool.py
GTGL = 'GL'  # GT else: GL
CHK_SZ = 10000
SUBCHUNK = 1000
WD = os.path.join(DATA_PATH)

PATH_GT_FILES = os.path.join(DATA_PATH, 'gt', 'stratified')
PATH_GL_FILES = os.path.join(DATA_PATH, 'gl')

SRCFILE = 'ALL.chr20.snps.gt.vcf.gz'
PATH_OUT = ['ALL.chr20.pooled.snps.{}.chunk{}.vcf'.format(GTGL.lower(), CHK_SZ),
            'ALL.chr20.missing.snps.{}.chunk{}.vcf'.format(GTGL.lower(), CHK_SZ)]
MSS = [False, True]
POOL = [True, False]
CHKFILE = 'ALL.chr20.snps.gt.chunk{}.vcf.gz'.format(CHK_SZ)

# unknown_gl = [1/3, 1/3, 1/3]
# unknown_gl = [0.2, 0.4, 0.4]
# unknown_gl = [0.02, 0.49, 0.49]
# unknown_gl = [0.02, 0.59, 0.39]
# unknown_gl = [0.02, 0.39, 0.59]
# unknown_gl = [0.02, 0.501, 0.479]
# unknown_gl = [0.02, 0.51, 0.47]
unknown_gl = 'adaptative'



### beagle.py
BEAGLE_JAR = os.path.join(SCRIPTS_PATH, 'beagle.11Mar19.69c.jar')
CFGT_JAR = os.path.join(SCRIPTS_PATH, 'conform-gt.jar')

RAW = {'vcf':'ALL.chr20.snps.{}.chunk{}.vcf'.format('gt', CHK_SZ),
       'gz':'ALL.chr20.snps.{}.chunk{}.vcf.gz'.format('gt', CHK_SZ),
       'ref': 'REF.chr20.snps.{}.chunk{}.vcf.gz'.format('gt', CHK_SZ),
       'imp': 'IMP.chr20.snps.{}.chunk{}.vcf.gz'.format('gt', CHK_SZ),
       'b1r':'REF.chr20.beagle1.chunk{}'.format(CHK_SZ),
       'b1i':'IMP.chr20.beagle1.chunk{}'.format(CHK_SZ)}

POOLED = {'vcf':'ALL.chr20.pooled.snps.{}.chunk{}.vcf'.format(GTGL.lower(), CHK_SZ),
          'gz':'ALL.chr20.pooled.snps.{}.chunk{}.vcf.gz'.format(GTGL.lower(), CHK_SZ),
          'imp': 'IMP.chr20.pooled.snps.{}.chunk{}.vcf.gz'.format(GTGL.lower(), CHK_SZ),
          'ref': 'REF.chr20.pooled.snps.{}.chunk{}.vcf.gz'.format(GTGL.lower(), CHK_SZ), # for MAF/AAF comparisons
          'b1':'IMP.chr20.pooled.beagle1.chunk{}'.format(CHK_SZ),
          'b2':'IMP.chr20.pooled.beagle2.{}.chunk{}'.format(GTGL.lower(), CHK_SZ),
          'corr':'IMP.chr20.pooled.beagle2.{}.chunk{}.corr'.format(GTGL.lower(), CHK_SZ),
          'cfgt': 'IMP.chr20.pooled.cfgt.chunk{}'.format(CHK_SZ),
          'gtonly': 'IMP.chr20.pooled.imputed.gt.chunk{}'.format(CHK_SZ)}

MISSING = {'vcf':'ALL.chr20.missing.snps.{}.chunk{}.vcf'.format(GTGL.lower(), CHK_SZ),
           'gz':'ALL.chr20.missing.snps.{}.chunk{}.vcf.gz'.format(GTGL.lower(), CHK_SZ),
           'imp': 'IMP.chr20.missing.snps.{}.chunk{}.vcf.gz'.format(GTGL.lower(), CHK_SZ),
           'b1':'IMP.chr20.missing.beagle1.chunk{}'.format(CHK_SZ),
           'b2':'IMP.chr20.missing.beagle2.{}.chunk{}'.format(GTGL.lower(), CHK_SZ),
           'corr': 'IMP.chr20.missing.beagle2.{}.chunk{}.corr'.format(GTGL.lower(), CHK_SZ),
           'cfgt': 'IMP.chr20.missing.cfgt.chunk{}'.format(CHK_SZ),
           'gtonly': 'IMP.chr20.missing.imputed.gt.chunk{}'.format(CHK_SZ)}

# To check: related individuals are removed from the file
try:
    nb_samples = len(VCF(os.path.join(WD,
                                      'gt',
                                      'ALL.chr20.snps.gt.chunk{}.vcf.gz'.format(CHK_SZ))).samples)
except IOError:
    nb_samples = np.nan
pools_size = 16
# number of pools for the target set IMP
pools_imp = ppcm(nb_samples, pools_size) // (10 * pools_size)
# number of samples for the target/reference sets
NB_IMP = pools_imp * pools_size
NB_REF = nb_samples - NB_IMP


### subchunking experiment
SUBSET = False
BIN_AAF = False
INTER = np.arange(0, 1, 0.1)


### PRINTING
print('*'.ljust(80, '*'))
print('Parameters configured:')
print('working directory: ', WD)
print('number of markers to extract and process: ', CHK_SZ)
print('number of pools: ', nb_samples // pools_size)
print('subsetting the main data set at 10%: ', SUBSET)
print('kind of genotype data: ', GTGL)
print('Path to GT files: ', PATH_GT_FILES)
print('Path to GL files: ', PATH_GL_FILES)
print('binarize MAFs: ', BIN_AAF)
if GTGL == 'GL':
    print('values for missing GLs: ', unknown_gl)
print('*'.ljust(80, '*'))
print('\r\n')