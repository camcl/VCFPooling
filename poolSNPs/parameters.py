#!/usr/bin/env python
"""
Parameters for running pooling + decoding on a dataset:
- path to file to process
- make new directory if needed: today's date, split data and pictures?
- files names/extension for outputs
- iteration number missing/pooled/missing.pooled
- chunks size if needed
"""

from cyvcf2 import  VCF

global CHK_SZ, SOURCE, KIND, MSS, POOL, WD
global RAW, POOLED, MISSING, MISS_POOL

### pool.py
WD = '/home/camilleclouard/PycharmProjects/1000Genomes/data/tests-beagle'
CHK_SZ = 1000
PATH_IN = 'ALL.chr20.snps.gt.vcf.gz'
PATH_OUT = ['ALL.chr20.pooled.snps.gt.chunk{}.vcf'.format(str(CHK_SZ)),
            'ALL.chr20.pooled.missing.snps.gt.chunk{}.vcf'.format(str(CHK_SZ)),
            'ALL.chr20.missing.snps.gt.chunk{}.vcf'.format(str(CHK_SZ))]
KIND = 'gt'
MSS = [False, True, True]
POOL = [True, True, False]
SOURCE = 'ALL.chr20.snps.gt.chunk{}.vcf.gz'.format(str(CHK_SZ))


### beagle.py
RAW = {'vcf':'ALL.chr20.snps.gt.chunk{}.vcf'.format(str(CHK_SZ)),
       'gz':'ALL.chr20.snps.gt.chunk{}.vcf.gz'.format(str(CHK_SZ)),
       'ref': 'REF.chr20.snps.gt.chunk{}.vcf.gz'.format(str(CHK_SZ)),
        'imp': 'IMP.chr20.snps.gt.chunk{}.vcf.gz'.format(str(CHK_SZ)),
       'b1r':'REF.chr20.beagle1',
       'b1i':'IMP.chr20.beagle1'}

POOLED = {'vcf':'ALL.chr20.pooled.snps.gt.chunk{}.vcf'.format(str(CHK_SZ)),
       'gz':'ALL.chr20.pooled.snps.gt.chunk{}.vcf.gz'.format(str(CHK_SZ)),
       'imp': 'IMP.chr20.pooled.snps.gt.chunk{}.vcf.gz'.format(str(CHK_SZ)),
       'b1':'IMP.chr20.pooled.beagle1',
       'b2':'IMP.chr20.pooled.beagle2',
          'corr':'IMP.chr20.pooled.beagle2.corr',
        'cfgt': 'IMP.chr20.pooled.cfgt'}

MISSING = {'vcf':'ALL.chr20.missing.snps.gt.chunk{}.vcf'.format(str(CHK_SZ)),
       'gz':'ALL.chr20.missing.snps.gt.chunk{}.vcf.gz'.format(str(CHK_SZ)),
       'imp': 'IMP.chr20.missing.snps.gt.chunk{}.vcf.gz'.format(str(CHK_SZ)),
       'b1':'IMP.chr20.missing.beagle1',
       'b2':'IMP.chr20.missing.beagle2',
           'corr': 'IMP.chr20.missing.beagle2.corr',
           'cfgt': 'IMP.chr20.missing.cfgt'}

MISS_POOL = {'vcf':'ALL.chr20.missing.snps.gt.chunk{}.vcf'.format(str(CHK_SZ)),
       'gz':'ALL.chr20.missing.pooled.snps.gt.chunk{}.vcf.gz'.format(str(CHK_SZ)),
       'imp': 'IMP.chr20.missing.pooled.snps.gt.chunk{}.vcf.gz'.format(str(CHK_SZ)),
       'b1':'IMP.chr20.missing.pooled.beagle1',
       'b2':'IMP.chr20.missing.pooled.beagle2',
             'corr': 'IMP.chr20.missing.pooled.beagle2.corr',
             'cfgt': 'IMP.chr20.missing.pooled.cfgt'}