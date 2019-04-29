#!/usr/bin/env python

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


### GENERAL TOOLS


def pgcd(a, b):
    """pgcd(a,b): calcul du 'Plus Grand Commun Diviseur' entre les 2 nombres entiers a et b"""
    while b != 0:
        r = a % b
        a, b = b, r
    return a


def ppcm(a, b):
    """ppcm(a,b): calcul du 'Plus Petit Commun Multiple' entre 2 nombres entiers a et b"""
    if (a == 0) or (b == 0):
        return 0
    else:
        return (a*b)//pgcd(a, b)


ext_vcf = '.vcf'
ext_vcfgz = '.vcf.gz'
ext_csv = '.csv'


### all
DATA_PATH = '/home/camille/PycharmProjects/1000Genomes/data'
PLOTS_PATH = '/home/camille/PycharmProjects/1000Genomes/plots'
SCRIPTS_PATH = '/home/camille/PycharmProjects/1000Genomes/scripts'


### pool.py
WD = os.path.join(DATA_PATH, 'tests-beagle')
CHK_SZ = 10000
SUBCHUNK = 1000
PATH_IN = 'ALL.chr20.snps.gt.vcf.gz'
PATH_OUT = ['ALL.chr20.pooled.snps.gt.chunk{}.vcf'.format(CHK_SZ),
            'ALL.chr20.missing.snps.gt.chunk{}.vcf'.format(CHK_SZ)]
#'ALL.chr20.pooled.missing.snps.gt.chunk{}.vcf'.format(CHK_SZ),

MSS = [False, True]
POOL = [True, False]
SOURCE = 'ALL.chr20.snps.gt.chunk{}.vcf.gz'.format(CHK_SZ)


### beagle.py
GTGL = 'GL'  # else: GL
BEAGLE_JAR = os.path.join(SCRIPTS_PATH, 'beagle.11Mar19.69c.jar')
CFGT_JAR = os.path.join(SCRIPTS_PATH, 'conform-gt.jar')

RAW = {'vcf':'ALL.chr20.snps.{}.chunk{}.vcf'.format(GTGL.lower(), CHK_SZ),
       'gz':'ALL.chr20.snps.{}.chunk{}.vcf.gz'.format(GTGL.lower(), CHK_SZ),
       'ref': 'REF.chr20.snps.{}.chunk{}.vcf.gz'.format(GTGL.lower(), CHK_SZ),
       'imp': 'IMP.chr20.snps.{}.chunk{}.vcf.gz'.format(GTGL.lower(), CHK_SZ),
       'b1r':'REF.chr20.beagle1.chunk{}'.format(CHK_SZ),
       'b1i':'IMP.chr20.beagle1.chunk{}'.format(CHK_SZ)}

POOLED = {'vcf':'ALL.chr20.pooled.snps.{}.chunk{}.vcf'.format(GTGL.lower(), CHK_SZ),
          'gz':'ALL.chr20.pooled.snps.{}.chunk{}.vcf.gz'.format(GTGL.lower(), CHK_SZ),
          'imp': 'IMP.chr20.pooled.snps.{}.chunk{}.vcf.gz'.format(GTGL.lower(), CHK_SZ),
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

# unknown_gl = [1/3, 1/3, 1/3]
# unknown_gl = [0.2, 0.4, 0.4]
unknown_gl = [0.02, 0.49, 0.49]
# unknown_gl = [0.02, 0.59, 0.39]
# unknown_gl = [0.02, 0.39, 0.59]
# unknown_gl = [0.02, 0.501, 0.479]
# unknown_gl = [0.02, 0.51, 0.47]
# unknown_gl = 'adaptative'




### subchunking experiment
SUBSET = False
BIN_MAF = False
INTER = np.arange(0, 1, 0.1)


### PRINTING
print('*'.ljust(80, '*'))
print('Parameters configured:')
print('working directory: ', WD)
print('number of markers to extract and process: ', CHK_SZ)
print('number of pools: ', nb_samples // pools_size)
print('subsetting the main data set at 10%: ', SUBSET)
print('kind of genotype data: ', GTGL)
print('binarize MAFs: ', BIN_MAF)
if GTGL == 'GL':
    print('values for missing GLs: ', unknown_gl)
print('*'.ljust(80, '*'))
print('\r\n')