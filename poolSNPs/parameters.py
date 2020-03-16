#!/usr/bin/env python3

#TODO: Rename files with core + extension (.vcf.gz, .vcf, .csv, ...)
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
# make the home path machine independent (no reconfiguring needed when running from server)
HOME = str(Path.home())
# TODO: replace all occurrences of /home/camille by prm.HOME
ROOT = os.path.join(HOME, '1000Genomes')
SNIC_PROJ = '/crex/proj/snic2019-8-216'

### all
DATA_PATH = os.path.join(ROOT, 'data')
TMP_DATA_PATH = os.path.join(ROOT, 'tmp_data')
PLOTS_PATH = os.path.join(ROOT, 'plots')
SCRIPTS_PATH = os.path.join(ROOT, 'scripts')


### pool.py
GTGL = 'GL'  # GT else: GL
CHK_SZ = None
if CHK_SZ is not None:
    chk_name = '.chunk{}'.format(CHK_SZ)
else:
    chk_name = ''
SUBCHUNK = 1000
WD = os.path.join(DATA_PATH)

PATH_GT_FILES = os.path.join(DATA_PATH, 'gt', 'stratified')
PATH_GL_FILES = os.path.join(DATA_PATH, 'gl')

SRCFILE = 'ALL.chr20.snps.gt.vcf.gz'
# PATH_OUT = ['ALL.chr20.pooled.snps.{}.chunk{}.vcf'.format(GTGL.lower(), CHK_SZ),
#             'ALL.chr20.missing.snps.{}.chunk{}.vcf'.format(GTGL.lower(), CHK_SZ)]
PATH_OUT = {'pooled': 'ALL.chr20.pooled.snps.{}.chunk{}.vcf'.format(GTGL.lower(), chk_name),
            'missing': 'ALL.chr20.missing.snps.{}.chunk{}.vcf'.format(GTGL.lower(), chk_name)}
MSS = [False, True]
POOL = [True, False]
CHKFILE = 'ALL.chr20.snps.gt.chunk{}.vcf.gz'.format(chk_name)

# unknown_gl = [1/3, 1/3, 1/3]
# unknown_gl = [0.2, 0.4, 0.4]
# unknown_gl = [0.02, 0.49, 0.49]
# unknown_gl = [0.02, 0.59, 0.39]
# unknown_gl = [0.02, 0.39, 0.59]
# unknown_gl = [0.02, 0.501, 0.479]
# unknown_gl = [0.02, 0.51, 0.47]
unknown_gl = 'adaptive'


### beagle.py
BEAGLE_JAR = os.path.join(SCRIPTS_PATH, 'beagle.11Mar19.69c.jar')
CFGT_JAR = os.path.join(SCRIPTS_PATH, 'conform-gt.jar')
MAP_FILE = os.path.join(DATA_PATH, 'plink.GRCh37.map', 'plink.chr20.GRCh37.map')

#TODO: rename beagle2.corr to imputed.gtdsgp
#TODO: rename IMP to STU and REF to PAN

RAW = {'vcf':'ALL.chr20.snps.{}{}.vcf'.format('gt', chk_name),
       'gz':'ALL.chr20.snps.{}{}.vcf.gz'.format('gt', chk_name),
       'ref': 'REF.chr20.snps.{}{}.vcf.gz'.format('gt', chk_name),
       'imp': 'IMP.chr20.snps.{}{}.vcf.gz'.format('gt', chk_name),
       'b1r':'REF.chr20.beagle1{}'.format(chk_name),
       'b1i':'IMP.chr20.beagle1{}'.format(chk_name)}

POOLED = {'vcf':'ALL.chr20.pooled.snps.{}{}.vcf'.format(GTGL.lower(), chk_name),
          'gz':'ALL.chr20.pooled.snps.{}{}.vcf.gz'.format(GTGL.lower(), chk_name),
          'imp': 'IMP.chr20.pooled.snps.{}{}.vcf.gz'.format(GTGL.lower(), chk_name),
          'ref': 'REF.chr20.pooled.snps.{}{}.vcf.gz'.format(GTGL.lower(), chk_name),  # for MAF/AAF comparisons
          'b1':'IMP.chr20.pooled.beagle1{}'.format(chk_name),
          'b2':'IMP.chr20.pooled.beagle2.{}{}'.format(GTGL.lower(), chk_name),
          'corr':'IMP.chr20.pooled.beagle2.{}{}.corr'.format(GTGL.lower(), chk_name),
          'cfgt': 'IMP.chr20.pooled.cfgt{}'.format(chk_name),
          'gtonly': 'IMP.chr20.pooled.imputed.gt{}'.format(chk_name)}

MISSING = {'vcf':'ALL.chr20.missing.snps.{}{}.vcf'.format(GTGL.lower(), chk_name),
           'gz':'ALL.chr20.missing.snps.{}{}.vcf.gz'.format(GTGL.lower(), chk_name),
           'imp': 'IMP.chr20.missing.snps.{}{}.vcf.gz'.format(GTGL.lower(), chk_name),
           'b1':'IMP.chr20.missing.beagle1{}'.format(chk_name),
           'b2':'IMP.chr20.missing.beagle2.{}{}'.format(GTGL.lower(), chk_name),
           'corr': 'IMP.chr20.missing.beagle2.{}{}.corr'.format(GTGL.lower(), chk_name),
           'cfgt': 'IMP.chr20.missing.cfgt{}'.format(chk_name),
           'gtonly': 'IMP.chr20.missing.imputed.gt{}'.format(chk_name)}

# To check: related individuals are removed from the file
try:
    nb_samples = len(VCF(os.path.join(WD,
                                      'gt',
                                      'ALL.chr20.snps.gt{}.vcf.gz'.format(chk_name))).samples)
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
def info():
    print('\n'.ljust(80, '*'))
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
    print('\r\n'.rjust(80, '*'))


def getter() -> dict:
    params = dict()
    params.update([('HOME', HOME), ('ROOT', ROOT), ('SNIC_PROJ', SNIC_PROJ)])
    params.update([('DATA_PATH', DATA_PATH), ('TMP_DATA_PATH', TMP_DATA_PATH),
                   ('SCRIPTS_PATH', SCRIPTS_PATH), ('PLOTS_PATH', PLOTS_PATH)])
    params.update([('WD', WD), ('PATH_GT_FILES', PATH_GT_FILES), ('PATH_GL_FILES', PATH_GL_FILES)])
    params.update([('SRCFILE', SRCFILE), ('GTGL', GTGL), ('unknown_gl', unknown_gl)
                   ('CHK_SZ', CHK_SZ), ('PATH_OUT', PATH_OUT), ('CHKFILE', CHKFILE)])
    params.update([('BEAGLE_JAR', BEAGLE_JAR), ('CFGT_JAR', CFGT_JAR), ('MAP_FILE', MAP_FILE)])
    params.update([('RAW', RAW), ('POOLED', POOLED), ('MISSING', MISSING)])
    params.update([('NB_IMP', NB_IMP), ('NB_REF', NB_REF)])
    params.update([('GLtype', GLtype), ('GTtype', GTtype)])
    return params
