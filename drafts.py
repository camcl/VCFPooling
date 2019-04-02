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
import pickle
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


vcf = VCF('IMP.chr20.pooled.beagle2.gt.chunk1000.corr.vcf.gz')

gl1000 = 'IMP.chr20.pooled.snps.gl.chunk1000.vcf.gz'
gt1000 = 'IMP.chr20.pooled.snps.gt.chunk1000.vcf.gz'
print(os.path.isfile(gl1000))

vcfgl = VCF(gl1000)
vcfgt = VCF(gt1000)

for vgl, vgt in zip(vcfgl('20:6304002'), vcfgt('20:6304002')):
    print(vgl)
    print(vgt.genotypes)
    print(alltls.bin_gl_converter(vgt))

alltls.file_likelihood_converter(gt1000, 'bug_gl1000.vcf')

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