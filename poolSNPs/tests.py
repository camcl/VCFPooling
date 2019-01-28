import os
import time
import warnings
import numpy as np
import pandas as pd
import math
from cyvcf2 import VCF, Writer
import itertools
import inspect

from scripts.poolSNPs import pool

wd = '/home/camilleclouard/PycharmProjects/1000Genomes/data/tests-beagle'
os.chdir(wd)

## Mtx test
M0 = {'in': np.array([[[0, 0, 1], [0, 0, 1], [0, 0, 1], [0, 0, 1]],
               [[0, 0, 1], [0, 0, 1], [0, 0, 1], [0, 0, 1]],
               [[0, 0, 1], [0, 0, 1], [0, 0, 1], [0, 0, 1]],
               [[0, 0, 1], [0, 0, 1], [0, 0, 1], [0, 0, 1]]]),
      'out': np.array([[[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]],
               [[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]],
               [[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]],
               [[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]]])}

M1 = {'in': np.array([[[1, 0, 1], [0, 0, 1], [0, 0, 1], [0, 0, 1]],
               [[0, 0, 1], [0, 0, 1], [0, 0, 1], [0, 0, 1]],
               [[0, 0, 1], [0, 0, 1], [0, 0, 1], [0, 0, 1]],
               [[0, 0, 1], [0, 0, 1], [0, 0, 1], [0, 0, 1]]]),
      'out': np.array([[[1, -1, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]],
               [[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]],
               [[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]],
               [[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]]])}

M2 = {'in': np.array([[[1, 1, 1], [0, 0, 1], [0, 0, 1], [0, 0, 1]],
               [[0, 0, 1], [0, 0, 1], [0, 0, 1], [0, 0, 1]],
               [[0, 0, 1], [0, 0, 1], [0, 0, 1], [0, 0, 1]],
               [[0, 0, 1], [0, 0, 1], [0, 0, 1], [0, 0, 1]]]),
      'out': np.array([[[1, -1, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]],
               [[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]],
               [[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]],
               [[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]]])}

M3 = {'in': np.array([[[1, 0, 1], [0, 0, 1], [0, 0, 1], [0, 0, 1]],
               [[0, 0, 1], [0, 0, 1], [1, 0, 1], [0, 0, 1]],
               [[0, 0, 1], [0, 0, 1], [0, 0, 1], [0, 0, 1]],
               [[0, 0, 1], [0, 0, 1], [0, 0, 1], [0, 0, 1]]]),
      'out': np.array([[[-1, -1, 0], [0, 0, 0], [-1, -1, 0], [0, 0, 0]],
               [[-1, -1, 0], [0, 0, 0], [-1, -1, 0], [0, 0, 0]],
               [[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]],
               [[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]]])}

M4 = {'in': np.array([[[1, 1, 1], [0, 0, 1], [0, 0, 1], [0, 0, 1]],
               [[0, 0, 1], [0, 0, 1], [1, 0, 1], [0, 0, 1]],
               [[0, 0, 1], [0, 0, 1], [0, 0, 1], [0, 0, 1]],
               [[0, 0, 1], [0, 0, 1], [0, 0, 1], [0, 0, 1]]]),
      'out': np.array([[[-1, -1, 0], [0, 0, 0], [-1, -1, 0], [0, 0, 0]],
               [[-1, -1, 0], [0, 0, 0], [-1, -1, 0], [0, 0, 0]],
               [[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]],
               [[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]]])}

M5 = {'in': np.array([[[1, 1, 1], [0, 0, 1], [0, 0, 1], [0, 0, 1]],
               [[1, 1, 1], [0, 0, 1], [0, 0, 1], [0, 0, 1]],
               [[1, 1, 1], [0, 0, 1], [0, 0, 1], [0, 0, 1]],
               [[1, 1, 1], [0, 0, 1], [0, 0, 1], [0, 0, 1]]]),
      'out': np.array([[[1, 1, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]],
               [[1, 1, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]],
               [[1, 1, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]],
               [[1, 1, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]]])}

M6 = {'in': np.array([[[1, 1, 1], [0, 0, 1], [0, 0, 1], [0, 0, 1]],
               [[1, 1, 1], [0, 0, 1], [1, 0, 1], [0, 0, 1]],
               [[1, 1, 1], [0, 0, 1], [0, 0, 1], [0, 0, 1]],
               [[1, 1, 1], [0, 0, 1], [0, 0, 1], [0, 0, 1]]]),
      'out': np.array([[[1, 1, 0], [0, 0, 0], [-1, -1, 0], [0, 0, 0]],
               [[1, 1, 0], [0, 0, 0], [-1, -1, 0], [0, 0, 0]],
               [[1, 1, 0], [0, 0, 0], [-1, -1, 0], [0, 0, 0]],
               [[1, 1, 0], [0, 0, 0], [-1, -1, 0], [0, 0, 0]]])}

M7 = {'in': np.array([[[1, 1, 1], [1, 1, 1], [1, 1, 1], [1, 1, 1]],
               [[1, 1, 1], [0, 0, 1], [1, 0, 1], [0, 0, 1]],
               [[1, 1, 1], [0, 0, 1], [0, 0, 1], [0, 0, 1]],
               [[1, 1, 1], [0, 0, 1], [0, 0, 1], [0, 0, 1]]]),
      'out': np.array([[[1, 1, 0], [1, 1, 0], [1, 1, 0], [1, 1, 0]],
               [[1, 1, 0], [-1, -1, 0], [-1, -1, 0], [-1, -1, 0]],
               [[1, 1, 0], [-1, -1, 0], [-1, -1, 0], [-1, -1, 0]],
               [[1, 1, 0], [-1, -1, 0], [-1, -1, 0], [-1, -1, 0]]])}

M8 = {'in': np.array([[[1, 1, 1], [1, 1, 1], [1, 1, 1], [1, 1, 1]],
               [[1, 1, 1], [1, 1, 1], [1, 1, 1], [1, 1, 1]],
               [[1, 1, 1], [1, 1, 1], [1, 1, 1], [1, 1, 1]],
               [[1, 1, 1], [1, 1, 1], [1, 1, 1], [1, 1, 1]]]),
      'out': np.array([[[1, 1, 0], [1, 1, 0], [1, 1, 0], [1, 1, 0]],
               [[1, 1, 0], [1, 1, 0], [1, 1, 0], [1, 1, 0]],
               [[1, 1, 0], [1, 1, 0], [1, 1, 0], [1, 1, 0]],
               [[1, 1, 0], [1, 1, 0], [1, 1, 0], [1, 1, 0]]])}

M9 = {'in': np.array([[[-1, -1, 1], [-1, -1, 1], [-1, -1, 1], [-1, -1, 1]],
               [[-1, -1, 1], [-1, -1, 1], [-1, -1, 1], [-1, -1, 1]],
               [[-1, -1, 1], [-1, -1, 1], [-1, -1, 1], [-1, -1, 1]],
               [[-1, -1, 1], [-1, -1, 1], [-1, -1, 1], [-1, -1, 1]]]),
      'out': np.array([[[-1, -1, 0], [-1, -1, 0], [-1, -1, 0], [-1, -1, 0]],
               [[-1, -1, 0], [-1, -1, 0], [-1, -1, 0], [-1, -1, 0]],
               [[-1, -1, 0], [-1, -1, 0], [-1, -1, 0], [-1, -1, 0]],
               [[-1, -1, 0], [-1, -1, 0], [-1, -1, 0], [-1, -1, 0]]])}


def pool_genotypes(d, call):
    """
    Computes genotypes of the different pools.
    :return: array of {0,1} GLs for RR|RA|AA for each pool
    """
    pools_size = 8
    scores = np.apply_along_axis(np.sum(call[:-1], axis=1))

    if np.any(call, -1):
        p = np.asarray([-1, -1, 0]*(len(call)//2))
    else:
        pooled_gt = np.dot(d, np.transpose(scores))
        bin_gt = np.vectorize(lambda x: 0 if x == 0 else (2 if x == 2*pools_size else 1))
        p = bin_gt(pooled_gt)
    return p # list of gt for the 8 pools from design matrix


def decoded_genotypes(pooled_samples, pools_list, pooled, subset, drop=False):
    """
    Recoomputes genotypes of samples with/without pooling/missing data
    :param pooled_samples: Variant.genotypes
    :return:
    """
     # = self.pools_list()
     # = self.pool_genotypes()
    nb_alt = np.sum(
        np.apply_along_axis(
            np.isin(pooled, 1)
        )
    )
    for s in subset: #16
        p = np.where(np.isin(np.asarray(pools_list), s))[0]
        out = pool.decoding_rules(pooled[:, p], nb_alt)
        pooled_samples[s] = out
        if drop == True:
            # missing data simulation
            np.put(pooled_samples[s],
                   [0, 1],
                   [-1, -1]
                   )

    return pooled_samples


def _pool_decode(M):
    Z = np.zeros_like(M['in']).reshape((1, 16, 3))
    oo = pool.SNPsPool()
    d = oo.design_matrix()
    f = M['in'].reshape((1, 16, 3))

    pools = np.zeros((1, 8, 3), dtype=int)
    for i in range(d.shape[0]):
        cache = d[i,:]
        poolval = f[:, np.argwhere(cache == 1)].squeeze()
        pools[:, i, :] = oo.pooling_rules(poolval)

    nb_alt = np.sum(
        np.apply_along_axis(
            lambda x: 1 in x,
            -1,
            pools[:, :, :-1]
        )
    )
    for i,s in enumerate(f[0]):  # 16
        p_p = np.argwhere(d[:,i] == 1).flatten()
        Z[:, i, :] = oo.decoding_rules(pools[:, p_p], nb_alt)

    return Z.reshape((4, 4, 3))

def test_pool_decode():
    # for i, m in enumerate([M0, M1, M2, M3, M4, M5, M6, M7, M8, M9]):
    for i, m in enumerate([M1]):
        z = _pool_decode(m)
        t = np.all(m['out'] == z)
        print(str(i) + ' --> ', t)
        if t == False:
            print(z)
        print('')

test_pool_decode()