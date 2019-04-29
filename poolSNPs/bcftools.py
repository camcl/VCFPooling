import subprocess
import os
import numpy as np

from scripts.poolSNPs import parameters as prm
from scripts.poolSNPs import patterns as pat
from scripts.poolSNPs.alleles import alleles_tools as alltls
from persotools.debugging import *
from persotools.files import *

"""
Bash commands for bcftools manipulations written as Python-functions.
Increase readability.
"""


def bgzip(f_vcf: str, f_gz: str, wd: str) -> None:
    """
    Bgzip a vcf file into a vcf.gz and checks if creation succeeded.
    :param f_vcf: input vcf file name
    :param f_gz: output vcf.gz file name
    :param wd: path to working directory
    :return: None
    """
    cmd = ' '.join(['bcftools',
                    'view -Oz -o',
                    f_gz,
                    f_vcf
                    ])
    subprocess.run(cmd, shell=True, cwd=wd)
    if os.path.exists(f_gz):
        print('File created: ', f_gz)


def sort(f_gz: str, wd: str) -> None:
    """
    Sorts a VCF-file per increasing marher position on chromosom.
    :param f_gz: input file name (vcf or vcf.gz)
    :param wd: path to working directory
    :return: None
    """
    cmd = ' '.join(['bcftools',
                    'sort -Oz -o',
                    f_gz,
                    f_gz
                    ])
    subprocess.run(cmd, shell=True, cwd=wd)


def index(f_gz: str, wd: str) -> None:
    """
    Index a VCF-file and replace old .csi file if it exixts.
    :param f_gz: input file name (vcf or vcf.gz)
    :param wd: path to working directory
    :return: None
    """
    cmd = ' '.join(['bcftools',
                    'index -f',
                    f_gz
                    ])
    subprocess.run(cmd, shell=True, cwd=wd)


def sampling(f_gz: str, f_out: str, f_samp: str, wd: str) -> None:
    """
    Creates a VCF-file (vcf.gz) for a chosen set of samples.
    Samples should be stored in a text-file, one sample name per line.
    Checks if sampled file was created.
    :param f_gz: input VCF-file name
    :param f_out: output VCF-file name (vcf.gz)
    :param f_samp: samples file name (text file)
    :param wd: path to working directory
    :return: None
    """
    cmd = ' '.join(['bcftools',
                    'view -Oz -o',
                    f_out,
                    '-S {}'.format(f_samp),
                    f_gz
                    ])
    subprocess.run(cmd, shell=True, cwd=wd)
    if os.path.exists(f_out):
        print('File created: ', f_out)


def stratified_maf_sampling(f_gz: str, wd: str) -> None:
    """
    Samples markers from a VCF-file by selecting an equal number of
    markers for each 0.1-MAF interval. Total number of markers depends
    on the global chunk size parameter.
    :param f_gz: input VCF-file with all markers
    :param wd: path to working directory
    :return: None
    """
    subprocess.run(' '.join(['bcftools',
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
        delete_file('ALL.chr20.snps.gt.chunk{}.strat.vcf'.format(prm.CHK_SZ))

