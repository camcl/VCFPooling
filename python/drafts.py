import os
import time
import warnings
import numpy as np
import math
import matplotlib.pyplot as plt
import itertools
from scripts.VCFPooling.poolSNPs.alleles import alleles_tools as alltls
from typing import *

from cyvcf2 import VCF
from pdbio.vcfdataframe import VcfDataFrame
import pysam

warnings.simplefilter('ignore')

start = time.time()
stop = time.time()
#print('Processing duration [s]: ', stop-start)


def grouper_it(n, iterable):
    it = iter(iterable)
    while True:
        chunk_it = itertools.islice(it, n)
        try:
            first_el = next(chunk_it)
        except StopIteration:
            return
        yield itertools.chain((first_el,), chunk_it)


os.chdir('/home/camille/1000Genomes/data/gl/gl_adaptive/all_snps_all_samples')
print(os.getcwd())


gl1000 = 'IMP.chr20.pooled.snps.gl.chunk1000.vcf.gz'
gt1000 = 'IMP.chr20.pooled.snps.gt.chunk1000.vcf.gz'
gl10000 = 'IMP.chr20.pooled.snps.gl.chunk1000.vcf.gz'
gt10000 = 'IMP.chr20.pooled.snps.gt.chunk1000.vcf.gz'
gtgl10000 = 'IMP.chr20.pooled.beagle2.gl.chunk10000.vcf.gz'


def cyvcf2_iter():
    return VCF(gtgl10000)


def pdbio_df():
    return VcfDataFrame(path=gtgl10000)

def pysam_iter(fpath):
    return pysam.VariantFile(fpath)


def f1(x):
    return 0.5 + math.exp(-x) * (x - 0.5)**3


def f2(x):
    return 0.5 + (math.exp(-x-0.5) - 0.5)**3


def my_tuple(a: int, b: int, c:int) -> Tuple[int]:
    return tuple((a, b, c))


class Shadok(object):
    def __init__(self, name):
        self.name = name

    def assign_task(self, state):
        self.__setattr__('state', state)
        print('Shadok {} is now {}'.format(self.name, self.__getattribute__('state')))


if __name__ == '__main__':
    t47_5 = my_tuple(4, 7, -5)
    print(my_tuple.__annotations__)

    s1 = Shadok('shadok1')
    print(s1.__dict__)
    s1.assign_task('pumping')
    print(s1.__dict__)
    print(s1.name)
    print(s1.state)

    """
        #CHROM       POS  ...                     HG01500                     HG01058
0       20     60826  ...                     0,-5,-5                     0,-5,-5
1       20     61044  ...                     0,-5,-5                     0,-5,-5
2       20     61271  ...                     0,-5,-5                     0,-5,-5
3       20     61279  ...                     0,-5,-5                     0,-5,-5
4       20     61409  ...                     0,-5,-5                     0,-5,-5
..     ...       ...  ...                         ...                         ...
995     20  62596412  ...   -0.6994,-0.27014,-0.57948                     0,-5,-5
996     20  62692060  ...    -0.60879,-0.2944,-0.6087    -0.60879,-0.2944,-0.6087
997     20  62778921  ...   -0.60879,-0.2944,-0.60879   -0.60879,-0.2944,-0.60879
998     20  62846376  ...  -0.62567,-0.29106,-0.59926  -0.62567,-0.29106,-0.59926
999     20  62855194  ...  -0.62567,-0.29106,-0.59926  -0.62567,-0.29106,-0.59926
[1000 rows x 249 columns]
    #CHROM       POS  ...                     HG01500                     HG01058
0       20     60826  ...                     0,-5,-5                     0,-5,-5
1       20     61044  ...                     0,-5,-5                     0,-5,-5
2       20     61271  ...                     0,-5,-5                     0,-5,-5
3       20     61279  ...                     0,-5,-5                     0,-5,-5
4       20     61409  ...                     0,-5,-5                     0,-5,-5
..     ...       ...  ...                         ...                         ...
995     20  62596412  ...   -0.6994,-0.27014,-0.57948                     0,-5,-5
996     20  62692060  ...    -0.60879,-0.2944,-0.6087    -0.60879,-0.2944,-0.6087
997     20  62778921  ...   -0.60879,-0.2944,-0.60879   -0.60879,-0.2944,-0.60879
998     20  62846376  ...  -0.62567,-0.29106,-0.59926  -0.62567,-0.29106,-0.59926
999     20  62855194  ...  -0.62567,-0.29106,-0.59926  -0.62567,-0.29106,-0.59926
[1000 rows x 249 columns]
    """

    myvcf = pysam_iter(gtgl10000)
    for v in myvcf:
        var = v
        break

    print(var.samples.items())
    for k, g in var.samples.items():
        break
        print(k, g['GP'])

    mygt = pysam_iter('/home/camille/1000Genomes/data/gt/' + gt10000)
    print(list(mygt.header.samples))
    missing = np.vectorize(lambda x: np.nan if x is None else x)
    for v in mygt:
        break
        genotypes = np.array([g['GT'] for g in v.samples.values()])
        tri = missing(genotypes.astype(float)).sum(axis=-1)
        print(np.nan_to_num(tri, nan=-1))
