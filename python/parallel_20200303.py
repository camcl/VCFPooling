import sys, os
import numpy as np
import pandas as pd
import cyvcf2
from itertools import repeat
import shutil
import multiprocessing as mp
import argparse
import time

home_dir = os.path.expanduser("~")
proj_dir = os.path.join(home_dir, '1000Genomes')
sys.path.insert(0, proj_dir)

from src.VCFPooling.poolSNPs import parameters as prm
from src.VCFPooling.poolSNPs import chunkvcf as chkvcf
from src.VCFPooling.python.archived import pool
from src.VCFPooling.poolSNPs import pybcf

from persotools.files import delete_file

'''
Parallelized file processing
For VCF-file bigger than some dozen of thousands of varaiants (e.g. 1 million SNPs).

Steps:
* Read main vcf and write chunks: bash script based on bcftools
* Simulate pooling on chunks
* Merge pooled chunks back to read-for-use compressed VCF file

Usage: 
$ python3 parallel_20200303.py ~/1000Genomes/data/gt/ALL.chr20.snps.gt.chunk10000.vcf.gz ~/1000Genomes/data/gl/gl_adaptive/parallel_pooling/ALL.chr20.pooled.snps.gl.chunk10000.vcf.gz 2
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
    tmp_path = prm.TMP_DATA_PATH

try:
    assert os.path.exists(prm.SNIC_PROJ)
    snic_proj = prm.SNIC_PROJ
except AssertionError:
    snic_proj = prm.TMP_DATA_PATH
print('SNIC PROJ: {}'.format(snic_proj))

print('\n'.ljust(80, '*'))
print('Number of cpu to be used = {}'.format(nb_cores))
print('Input file = {}'.format(os.path.expanduser(argsin.pathin)))
print('Output file = {}'.format(os.path.expanduser(argsin.pathout)))
print('\n'.rjust(80, '*'))

### SPLIT MAIN VCF-FILE INTO PACKS
os.chdir(tmp_path)
fingz = os.path.expanduser(argsin.pathin)  # '/home/camille/...
basin = os.path.basename(fingz).rstrip('.gz')
basout = os.path.basename(argsin.pathout)

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

start = time.time()
prename = 'pack'

files0 = [f for f in os.listdir(os.path.join(snic_proj, 'parallel_vcfsplit'))]
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
#print([os.path.join(snic_proj, 'parallel_vcfsplit', f0) for f0 in files0])

files1 = ['pooled.{}'.format(f0).rstrip('.gz') for f0 in files0]
args1 = list(zip([os.path.join(snic_proj, 'parallel_vcfsplit', f0) for f0 in files0],
                 repeat(splits, len(indices)),
                 repeat(df2dict, len(indices)),
                 [os.path.join(tmp_path, f1) for f1 in files1]))
with mp.Pool(processes=nb_cores) as mpool:
    _ = all(mpool.starmap(chkvcf.cyvcfchunkhandler_process, args1))

print('\r\nTime elapsed --> ', time.time() - start)

files2 = ['{}.gz'.format(f1) for f1 in files1]
pybcf.concat(files2, basin, tmp_path)
pybcf.sort(basin, tmp_path)
pybcf.bgzip(basin, basout, tmp_path)
pybcf.index(basout, tmp_path)
delete_file(basin)
shutil.copy(basout, os.path.expanduser(argsin.pathout))
print('\r\nTime elapsed --> ', time.time() - start)
