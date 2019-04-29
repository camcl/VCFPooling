import os
import time
import warnings
from scipy.stats import bernoulli as bn
import numpy as np
import pandas as pd
from cyvcf2 import VCF, Writer
import itertools
import math
import inspect
import matplotlib.pyplot as plt
from operator import *
import itertools
from scripts.poolSNPs.alleles import alleles_tools as alltls


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


os.chdir('./data/tests-beagle')
print(os.getcwd())


gt1000 = 'IMP.chr20.snps.gt.chunk1000.vcf.gz'
df_maf = alltls.get_maf(gt1000)
dcount = df_maf['maf'].apply(lambda x: np.digitize(x,
                                                   np.arange(0, 1, 0.1))
                                       /10).to_frame(name='bin_maf')
#print(dcount.groupby(by='bin_maf', axis=0)['bin_maf'].count())
# dcount.plot.hist(by='maf_bin')
# plt.show()


gl1000 = 'IMP.chr20.pooled.snps.gl.chunk1000.vcf.gz'
gt1000 = 'IMP.chr20.pooled.snps.gt.chunk1000.vcf.gz'


# 4 slots, 2 idv
nbSlots = 4
nbIdv = 2
gtpos = ['RA', 'AA']
prd = itertools.combinations_with_replacement(gtpos, r=nbIdv)
patterns = list()
for p in prd:
    tpl = p + tuple(['RR'] * (nbSlots - nbIdv))
    patterns.append(tpl)
pot = np.array(patterns).flatten()
# weighted GLs
wgl = np.array([np.count_nonzero(np.where(pot == let, 1, 0)) for let in ['RR', 'RA', 'AA']]) / len(pot)
print(wgl)


"""
_x = np.arange(0, 1, 0.02)


def f1(x):
    return (1 / (1 + math.exp(-20*(x - 0.5))))


def f2(x):
    return (x**0.25 / (1 + math.exp(-20*(x - 0.5))))


_y1 = list(map(f1, _x))
_y2 = list(map(f2, _x))


fig, ax = plt.subplots()
ax.plot(_x.tolist(), _y1, c='g')
ax.plot(_x.tolist(), _y2, c='b')
ax.plot(_x.tolist(), _x.tolist(), c='k')
plt.show()
"""