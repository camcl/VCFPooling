import os
import argparse
import subprocess
import warnings
import allel
import numpy as np
import time
import multiprocessing as mp
import mmap

"""
README:
Script written to be used from the command line.

Requirements:
- bcftools
- Python 3.x: scikit-allel package (in the venv)

Usage:
/!\ prerequisite commands locally:
$ cd PycharmProjects/1000Genomes/scripts
$ source ../venv/bin/activate

Run for GLs:
$ python3 subsample.py ALL.chr20.phase3_GL_idv_ID.vcf ALL.chr20.phase3_GL_subset_ ALL.chr20.phase3_bc_union.20130502.biallelic_svm_snps_indels_ANNOTATED.vcf.gz 159
GTs:
$ python3 subsample.py ALL.chr20.snps.gt.txt ALL.chr20.snps.gt_subset_ ALL.chr20.snps.gt.vcf.gz 157
"""
start = time.time()

warnings.simplefilter('ignore')

os.chdir('../data')

# Parse the arguments for cmd usage
parser = argparse.ArgumentParser(description='''Randomly subsets the main VCF data into pools.
 Pools are written as new VCF files''')
parser.add_argument('idx',
                    action='store',
                    help='Name of the VCF file containing the individuals IDs')
parser.add_argument('sub',
                    action='store',
                    help='String prefix for naming the subset VCF files')
parser.add_argument('data',
                    action='store',
                    help='Name of the file containing the SNP calls/genetic data')
parser.add_argument('filesNb',
                    action='store',
                    help='Number of subset files to create')
args = parser.parse_args()

# Import the ID file
with open (args.idx, 'r') as f:
    idv_id = f.readlines()
idv_id = [idv.rstrip('\n') for idv in idv_id]

# Parameters
idv_nb = len(idv_id)
pool_size = 16
nb_pool, left = divmod(idv_nb, pool_size)
print('Number 16-sized pools to create, number of single samples reamaining: ',
      nb_pool, left)

# Create random pools
pools = np.empty((1,pool_size), dtype=np.str)
n = 0
while n < nb_pool:
    pools = np.vstack((pools, np.random.choice(idv_id, size=pool_size, replace=False)))
    if n == 0:
        pools = pools[1:,:]
    idv_id = np.delete(idv_id, np.argwhere(pools[n,:]))
    n+=1
singles = idv_id

# Create VCF pools files
try:
    GL = input('VCF file with GL? (T/F)').upper()
    while GL != 'F' and GL != 'T':
        GL = input('VCF file with GL? (T/F)').upper()
except ValueError:
    print('Enter the letter T or F')

# # Split pools set
# N = 2
# spl = np.array_split(pools[:4,:], N, axis=0) # rm pools slicing for processing the whole dataset
# m = max(len(a) for a in spl)
# idx = [0 + i*m for i in range(len(spl))]
#
#
# def subset(split, start_idx):
#     """
#
#     :param split:
#     :param start_idx:
#     :return:
#     """
#     for n, p in enumerate(range(len(split))):
#         n = start_idx + n
#         print(n)
#         print('Write VCF for pool {}'.format(str(n + 1).rjust(3, '0')).ljust(80, '.'))
#         if GL == 'F':
#             vcf_pl = ' '.join(['bcftools',
#                                'view -Oz',
#                                '-s ' + ','.join([i for i in split[p, :]]),
#                                '-o ' + args.sub + '{}.vcf.gz'.format(str(n + 1).rjust(3, '0')),
#                                args.data
#                                ])
#         else:
#             vcf_pl = ' '.join(['bcftools',
#                                '+fill-AN-AC',
#                                args.data,
#                                '|',
#                                'bcftools',
#                                'view -Oz',
#                                '-s ' + ','.join([i for i in split[p, :]]),
#                                '-o ' + args.sub + '{}.vcf.gz'.format(str(n + 1).rjust(3, '0'))
#                                ])
#
#         vcf_ix = ' '.join(['bcftools',
#                            'index',
#                            args.sub + '{}.vcf.gz'.format(str(n + 1).rjust(3, '0'))
#                            ])
#         subprocess.call([vcf_pl, vcf_ix], shell=True)
#     stop = time.time()
#    print('Elapsed time: {}'.format(str(stop - start)).rjust(80, '.'))


def subset2(n):
    """

    :param n: index of the subset to be processed
    :return:
    """
    print('Write VCF for pool {}'.format(str(n + 1).rjust(3, '0')).ljust(80, '.'))
    if GL == 'F':
        vcf_pl = ' '.join(['bcftools',
                           'view -Oz',
                           '-s ' + ','.join([i for i in pools[n, :]]),
                           '-o ' + args.sub + '{}.vcf.gz'.format(str(n + 1).rjust(3, '0')),
                           args.data
                           ])
    else:
        vcf_pl = ' '.join(['bcftools',
                           '+fill-AN-AC',
                           args.data,
                           '|',
                           'bcftools',
                           'view -Oz',
                           '-s ' + ','.join([i for i in pools[n, :]]),
                           '-o ' + args.sub + '{}.vcf.gz'.format(str(n + 1).rjust(3, '0'))
                           ])

    vcf_ix = ' '.join(['bcftools',
                       'index',
                       args.sub + '{}.vcf.gz'.format(str(n + 1).rjust(3, '0'))
                       ])
    subprocess.call([vcf_pl, vcf_ix], shell=True)


    stop = time.time()
    print('Elapsed time: {}'.format(str(stop - start)).rjust(80, '.'))

# Equivalents for running a single process
# subset(pools, 0)

# for n in range(nb_pool):
#     subset2(n)

def remaining():
    """
    Write VCF for the remaining individuals (not in any pool)
    :return:
    """
    if int(args.filesNb) == nb_pool + 1:
        print('Write VCF for the remaining individuals'.ljust(80, '.'))
        if GL == 'F':
            vcf_sg = ' '.join(['bcftools',
                               'view -Oz',
                               '-s ' + ','.join([i for i in singles]),
                               '-o ' + args.sub + '{}.vcf.gz'.format(str(nb_pool+1).rjust(3, '0')),
                               args.data
                          ])
        else:
            vcf_sg = ' '.join(['bcftools',
                               '+fill-AN-AC',
                               args.data,
                               '|',
                               'bcftools',
                               'view -Oz',
                               '-s ' + ','.join([i for i in singles]),
                               '-o ' + args.sub + '{}.vcf.gz'.format(str(nb_pool+1).rjust(3, '0'))
                          ])

        vcf_ix = ' '.join(['bcftools',
                           'index',
                           args.sub + '{}.vcf.gz'.format(str(nb_pool+1).rjust(3, '0'))
                      ])
        subprocess.call([vcf_sg, vcf_ix], shell=True)

        print('\r\nAll VCF subsets ready.')
        stop = time.time()
        print('Elapsed time: {}'.format(str(stop-start)).rjust(80, '.'))

### MULTIPROCESSING
# Test for 4 processes on 4 files
pll = mp.Pool(processes=4)
_ = pll.map(subset2, range(int(args.filesNb)-1)) # if filesNb = nb_pool, estimated time = 7960 sec * 40 = 318400 sec = 3 days 16:26:40
remaining()

# Test for " processes on 4 files
# param = spl, idx
# processes = [mp.Process(target=subset, args=(s, x)) for s,x in zip(spl, idx)]
#
# for pr in processes:
#     pr.start()
#
# for pr in processes:
#     pr.join()