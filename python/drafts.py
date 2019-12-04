import os
import time
import warnings
import numpy as np
import math
import matplotlib.pyplot as plt
import itertools
from scripts.VCFPooling.poolSNPs.alleles import alleles_tools as alltls
from typing import *


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

#print(dcount.groupby(by='bin_maf', axis=0)['bin_maf'].count())
# dcount.plot.hist(by='maf_bin')
# plt.show()


gl1000 = 'IMP.chr20.pooled.snps.gl.chunk1000.vcf.gz'
gt1000 = 'IMP.chr20.pooled.snps.gt.chunk1000.vcf.gz'



_x = np.arange(0, 1, 0.02)


def f1(x):
    return 0.5 + math.exp(-x) * (x - 0.5)**3


def f2(x):
    return 0.5 + (math.exp(-x-0.5) - 0.5)**3


_y1 = list(map(f1, _x))
_y2 = list(map(f2, _x))


fig, ax = plt.subplots()
ax.plot(_x.tolist(), _y1, c='g')
ax.plot(_x.tolist(), _y2, c='b')
ax.plot(_x.tolist(), _x.tolist(), c='k')


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

    print(alltls.SigmoidInterpolator)
    print(type(alltls.SigmoidInterpolator))


    s1 = Shadok('shadok1')
    print(s1.__dict__)
    s1.assign_task('pumping')
    print(s1.__dict__)
    print(s1.name)
    print(s1.state)
