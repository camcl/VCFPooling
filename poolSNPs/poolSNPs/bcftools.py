import subprocess
import os
import numpy as np

from scripts.poolSNPs import parameters as prm
from scripts.poolSNPs import patterns as pat
from scripts.poolSNPs.alleles import alleles_tools as alltls
from persotools.debugging import *
from persotools.files import *

'''
Commands for bcftools manipulations written as Python-functions.
Increase readability.
'''


def bgzip(f_vcf, f_gz, wd):
    cmd = ' '.join(['bcftools',
                    'view -Oz -o',
                    f_gz,
                    f_vcf
                    ])
    subprocess.run(cmd, shell=True, cwd=wd)
    if os.path.exists(f_gz):
        print('File created: ', f_gz)


def sort(f_gz, wd):
    cmd = ' '.join(['bcftools',
                    'sort -Oz -o',
                    f_gz,
                    f_gz
                    ])
    subprocess.run(cmd, shell=True, cwd=wd)


def index(f_gz, wd):
    cmd = ' '.join(['bcftools',
                    'index -f',
                    f_gz
                    ])
    subprocess.run(cmd, shell=True, cwd=wd)


def sampling(f_gz, f_out, f_samp, wd):
    cmd = ' '.join(['bcftools',
                    'view -Oz -o',
                    f_out,
                    '-S {}'.format(f_samp),
                    f_gz
                    ])
    subprocess.run(cmd, shell=True, cwd=wd)
    if os.path.exists(f_out):
        print('File created: ', f_out)


def stratified_maf_sampling(f_gz, wd):
    subprocess.run(''.join(['bcftools',
                            'view -h -Ov -o',
                            'ALL.chr20.snps.gt.chunk{}.strat.vcf'.format(prm.CHK_SZ),
                            'ALL.chr20.snps.gt.vcf.gz']),
                   shell=True, cwd=wd)
    bins = np.arange(0.0, 1.0, 0.1)
    for b in bins[:-1]:
        cmd1 = ' '.join(['bcftools',
                        'view -Oz -o',
                        f_gz.replace('.vcf.gz', '.maf_{}_{}.vcf.gz'.format(b, b+1)),
                        '-q {} -Q {}'.format(b, b+1),
                        f_gz
                        ])

        cmd2 = ' '.join(['bcftools view -H',
                         f_gz.replace('.vcf.gz', '.maf_{}_{}.vcf.gz'.format(b, b+1)),
                         '| sort -R | head -{}'.format(prm.CHK_SZ // len(bins)),
                         '| tee -a ALL.chr20.snps.gt.chunk{}.strat.vcf'.format(prm.CHK_SZ)
                         ])

        subprocess.run(cmd1, shell=True, cwd=wd)
        subprocess.run(cmd2, shell=True, cwd=wd)
        bgzip('ALL.chr20.snps.gt.chunk{}.strat.vcf'.format(prm.CHK_SZ),
              'ALL.chr20.snps.gt.chunk{}.strat.vcf.gz'.format(prm.CHK_SZ),
              wd)
        sort('ALL.chr20.snps.gt.chunk{}.strat.vcf.gz'.format(prm.CHK_SZ),
             wd)
        index('ALL.chr20.snps.gt.chunk{}.strat.vcf.gz'.format(prm.CHK_SZ),
             wd)

