import sys, os
from typing import *
import numpy as np
import pandas as pd
import cyvcf2
from itertools import starmap, repeat
import shutil
import multiprocessing as mp
import argparse
import time

home_dir = os.path.expanduser("~")
proj_dir = os.path.join(home_dir, '1000Genomes')
sys.path.insert(0, proj_dir)

from scripts.VCFPooling.poolSNPs import parameters as prm
from scripts.VCFPooling.poolSNPs import chunkvcf as chkvcf
from scripts.VCFPooling.poolSNPs import dataframe as vcfdf
from scripts.VCFPooling.poolSNPs import pool
from scripts.VCFPooling.poolSNPs import pybcf

from persotools.files import delete_file, mkdir
from persotools.struct import NamedDict

'''
Parallelized file processing

Steps:
* Read main vcf and write chunks
* Simulate pooling on chunks
* Merge pooled chunks back to read-for-use compressed VCF file

Usage: 
$ python3 parallel_20200222.py ~/1000Genomes/data/gt/ALL.chr20.snps.gt.chunk10000.vcf.gz ~/1000Genomes/data/gl/gl_adaptive/parallel_pooling/ALL.chr20.pooled.snps.gl.chunk10000.vcf.gz 2
'''

### COMMAND-LINE PARSING AND PARAMETERS
parser = argparse.ArgumentParser(description='Run parallelized pooling simulation'
                                             'on the whole set of samples')
parser.add_argument('pathin', metavar='in', type=str, help='File to pool', default=None)
parser.add_argument('pathout', metavar='out', type=str, help='Pooled file', default=None)
parser.add_argument('cores', metavar='nc', type=int, nargs='?', help='Number of cores to use', default=None)

argsin = parser.parse_args()

nb_cores = os.cpu_count() if argsin.cores is None else argsin.cores
try:
    tmp_pathh = os.environ['SNIC_TMP']
except KeyError:
    tmp_path = prm.TMP_DATA_PATH

print('\n'.ljust(80, '*'))
print('Number of cpu to be used = {}'.format(nb_cores))
print('Input file = {}'.format(argsin.pathin))
print('Output file = {}'.format(argsin.pathout))
print('\n'.rjust(80, '*'))

### SPLIT MAIN VCF-FILE INTO PACKS
os.chdir(tmp_path)
fingz = argsin.pathin  # '/home/camille/1000Genomes/data/gt/ALL.chr20.snps.gt.chunk10000.vcf.gz'
basin = os.path.basename(fingz).rstrip('.gz')
basout = os.path.basename(argsin.pathout)

cyvcfchunk = chkvcf.CyvcfVariantChunkGenerator(fingz, chunksize=1000)

df = pd.read_csv(os.path.join(prm.WD, 'adaptive_gls.csv'),
                 header=None,
                 names=['rowsrr', 'rowsra', 'rowsaa', 'colsrr', 'colsra', 'colsaa',
                        'n', 'm',
                        'rr', 'ra', 'aa']
                 )
df2dict = dict(((int(rwrr), int(rwra), int(rwaa), int(clrr), int(clra), int(claa),
                 int(n), int(m)),
                [rr, ra, aa]) for rwrr, rwra, rwaa, clrr, clra, claa,
                                  n, m,
                                  rr, ra, aa in df.itertuples(index=False, name=None))

raw = cyvcf2.VCF(os.path.join(prm.WD, 'gt', prm.SRCFILE))  # VCF iterator object
splits = pool.split_pools(raw.samples, 16, seed=123)  # list of lists
cyvcfpack = cyvcfchunk.chunkpacker()  # generator of generators

start = time.time()
prename = 'pack'
indices = np.arange(len([*cyvcfpack]))
files0 = [os.path.join(tmp_path,
                       '{}{}.{}'.format(prename, chki, basin)) for chki in indices]
args0 = list(zip(indices,
                 files0))
with mp.Pool(processes=os.cpu_count()) as mpool:
   _ = all(mpool.starmap(chkvcf.CyvcfVariantChunkGenerator(fingz, chunksize=1000).chunkwriter, args0))
print('\r\nTime elapsed --> ', time.time() - start)  # 12.9

files0 = os.listdir(tmp_path)
files1 = ['pooled.{}'.format(f0) for f0 in files0]
args1 = list(zip(files0,
                 repeat(splits, len(indices)),
                 repeat(df2dict, len(indices)),
                 files1))
with mp.Pool(processes=os.cpu_count()) as mpool:
     _ = all(mpool.starmap(chkvcf.cyvcfchunkhandler_process, args1))

files2 = ['{}.gz'.format(f1) for f1 in files1]
pybcf.concat(files2, basin, tmp_path)
pybcf.sort(basin, tmp_path)
pybcf.bgzip(basin, basout, tmp_path)
pybcf.index(basout, tmp_path)
delete_file(basin)
shutil.copy(basout, argsin.pathout)
print('\r\nTime elapsed --> ', time.time() - start)  # 19.45
