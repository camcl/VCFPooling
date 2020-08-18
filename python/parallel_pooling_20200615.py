import sys, os
import numpy as np
import pandas as pd
from itertools import starmap, repeat
import shutil
import subprocess
import multiprocessing as mp
import argparse
import timeit
import datetime

home_dir = os.path.expanduser("~")
bin_dir = os.path.join(home_dir, '1000Genomes')
proj_dir = os.path.join(home_dir, 'PoolImpHuman')
sys.path.insert(0, bin_dir)

from src.VCFPooling.poolSNPs import parameters as prm
from src.VCFPooling.poolSNPs import poolvcf
from src.VCFPooling.poolSNPs import pybcf

from persotools.files import delete_file, mkdir

'''
Parallelized file processing
For VCF-file bigger than some dozen of thousands of varaiants (e.g. 1 million SNPs).
Test here for 10,000 markers first.

Steps:
* Read main vcf and write chunks: bash script based on bcftools
* Simulate pooling on chunks
* Merge pooled chunks back to read-for-use compressed VCF file

Usage: 
$ python3 parallel_pooling_20200615.py ~/PoolImpHuman/data/20200615/ALL.chr20.snps.gt.chunk10000.vcf.gz ~/PoolImpHuman/data/20200615/ALL.chr20.pooled.snps.gl.chunk10000.vcf.gz 2
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
    tmp_path = os.environ['SNIC_TMP']
except KeyError:
    tmp_path = os.path.join(proj_dir, 'data/tmp')

try:
    assert os.path.exists(prm.SNIC_PROJ)
    snic_proj = '/crex/proj/snic2019-8-216/private'
except AssertionError:
    snic_proj = os.path.join(proj_dir, 'data')
print('SNIC PROJ: {}'.format(snic_proj))

print('\n'.ljust(80, '*'))
print('Number of cpu to be used = {}'.format(nb_cores))
print('Input file = {}'.format(os.path.expanduser(argsin.pathin)))
print('Output file = {}'.format(os.path.expanduser(argsin.pathout)))
print('\n'.rjust(80, '*'))

### SPLIT MAIN VCF-FILE INTO PACKS
today = datetime.datetime.today().date()
date_dir = today.strftime('%Y%m%d')
...
# Lines to run if directories and contents not created locally/remotely
# subprocess.run('mkdir {}'.format(date_dir), shell=True, cwd=snic_proj)
# subprocess.run('mkdir parallel_vcfsplit', shell=True, cwd=snic_proj)
# # bash bcfchunkpara.sh ...

### SIMULATE POOLING ON PACKED DATA CHUNKS
data_dir = prm.WD
os.chdir(tmp_path)
fingz = os.path.expanduser(argsin.pathin)  # /home/camille/...
basin = os.path.basename(fingz).rstrip('.gz')
basout = os.path.basename(argsin.pathout)

start = timeit.default_timer()
#prename = 'pack'

wd = os.path.join(snic_proj, 'parallel_vcfsplit')
files0 = [f for f in os.listdir(wd)]
try:
    files0.remove(basin + '.gz')
except ValueError:
    pass
try:
    files0.remove(basin + '.gz.csi')
except ValueError:
    pass
try:
    files0.remove(basin + '.gz.tbi')
except ValueError:
    pass

indices = np.arange(len(files0))

print('\r\n{} files found will be pooled'.format(len(files0)).ljust(80, '.'))
print([os.path.join(wd, f0) for f0 in files0])

files1 = ['pooled.{}'.format(f0).rstrip('.gz') for f0 in files0]
args1 = list(zip([os.path.join(wd, f0) for f0 in files0],
                 [os.path.join(tmp_path, f1) for f1 in files1],
                 repeat(data_dir, len(indices)),  # path to lookup table
                 repeat(tmp_path, len(indices))))  # directory for output

with mp.Pool(processes=nb_cores) as mpool:
    _ = all(mpool.starmap(poolvcf.pysam_pooler_gp, args1))

print('\r\nTime elapsed --> ', timeit.default_timer() - start)
# Time elapsed -->  588.0913308090021 (with 2 local cores)

### CONCATENATE RESULTS
files2 = ['{}.gz'.format(f1) for f1 in files1]
# files2 = ['/home/camille/PoolImpHuman/data/tmp/pooled.pack7.ALL.chr20.snps.gt.chunk10000.vcf.gz',
#           '/home/camille/PoolImpHuman/data/tmp/pooled.pack5.ALL.chr20.snps.gt.chunk10000.vcf.gz',
#           '/home/camille/PoolImpHuman/data/tmp/pooled.pack1.ALL.chr20.snps.gt.chunk10000.vcf.gz',
#           '/home/camille/PoolImpHuman/data/tmp/pooled.pack0.ALL.chr20.snps.gt.chunk10000.vcf.gz',
#           '/home/camille/PoolImpHuman/data/tmp/pooled.pack6.ALL.chr20.snps.gt.chunk10000.vcf.gz',
#           '/home/camille/PoolImpHuman/data/tmp/pooled.pack9.ALL.chr20.snps.gt.chunk10000.vcf.gz',
#           '/home/camille/PoolImpHuman/data/tmp/pooled.pack8.ALL.chr20.snps.gt.chunk10000.vcf.gz',
#           '/home/camille/PoolImpHuman/data/tmp/pooled.pack3.ALL.chr20.snps.gt.chunk10000.vcf.gz',
#           '/home/camille/PoolImpHuman/data/tmp/pooled.pack4.ALL.chr20.snps.gt.chunk10000.vcf.gz',
#           '/home/camille/PoolImpHuman/data/tmp/pooled.pack2.ALL.chr20.snps.gt.chunk10000.vcf.gz'
#           ]
pybcf.concat(files2, basin, tmp_path)
pybcf.sort(basin, tmp_path)
pybcf.bgzip(basin, basout, tmp_path)
pybcf.index(basout, tmp_path)
delete_file(basin)
shutil.copy(basout, os.path.expanduser(argsin.pathout))
print('\r\nTime elapsed --> ', timeit.default_timer() - start)
