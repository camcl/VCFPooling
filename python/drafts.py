import os
import time
import warnings
import numpy as np
import math
import matplotlib.pyplot as plt
import itertools
from typing import *

from cyvcf2 import VCF
import pysam
from scripts.VCFPooling.poolSNPs import _mylog
from scripts.VCFPooling.poolSNPs import _mypath

warnings.simplefilter('ignore')
_mylog.stdout(__file__)

from scripts.VCFPooling.poolSNPs.alleles import alleles_tools as alltls


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

gl1000 = 'IMP.chr20.pooled.snps.gl.chunk1000.vcf.gz'
gt1000 = 'IMP.chr20.pooled.snps.gt.chunk1000.vcf.gz'
gl10000 = 'IMP.chr20.pooled.snps.gl.chunk1000.vcf.gz'
gt10000 = 'IMP.chr20.pooled.snps.gt.chunk1000.vcf.gz'
gtgl10000 = 'IMP.chr20.pooled.beagle2.gl.chunk10000.vcf.gz'


def cyvcf2_iter():
    return VCF(gtgl10000)


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
    import timeit
    start_time = timeit.default_timer()

    os.chdir('/home/camille/1000Genomes/tmp_data')
    #print(bug)
    t47_5 = my_tuple(4, 7, -5)
    print(my_tuple.__annotations__)

    s1 = Shadok('shadok1')
    print(s1.__dict__)
    s1.assign_task('pumping')
    print(s1.__dict__)
    print(s1.name)
    print(s1.state)

    print('Process duration timeit: ', timeit.default_timer() - start_time)
    print('Process duration mix and mess: ', timeit.default_timer() - start)
    print('Process time.clock: ', time.clock())

    if False:
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
