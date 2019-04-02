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
from scripts.poolSNPs import parameters as prm
from scripts.poolSNPs.alleles import alleles_tools as alltls
from scripts.poolSNPs.alleles import alleles_plots as allplt


wd = '/home/camille/PycharmProjects/1000Genomes/data/tests-beagle'
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


def pool_genotypes(d, call):
    """
    Computes genotypes of the different pools.
    :return: array of {0,1} GLs for RR|RA|AA for each pool
    """
    scores = np.apply_along_axis(sum, axis=-1, arr=call[:, :, :-1]).reshape((1, 16, 1))

    if np.isin(call, -1).any():
        x = np.ones((1, 8, 1))
        y = np.asarray([-1, -1, 0])
        b = np.broadcast(x, y)
        p = np.empty(b.shape)
        p.flat = [u * v for (u, v) in b]
    else:
        pooled_gt = np.dot(d, np.transpose(scores)).reshape((1,8,1))
        pooled_gt = np.broadcast_to(pooled_gt, (1, 8, 3))
        pooler = lambda x: [0, 0, 0] if np.all(x == 0) else ([1, 1, 0] if np.all(x == 8) else [1, 0, 0])
        p = np.apply_along_axis(pooler, axis=-1, arr=pooled_gt)
    return p # list of gt for the 8 pools from design matrix


def decode_genotypes(d, call, drop=False):
    """
    Recoomputes genotypes of samples with/without pooling/missing data
    :param pooled_samples: Variant.genotypes
    :return:
    """
    pooled = pool_genotypes(d, call)
    scores = np.apply_along_axis(sum, axis=-1, arr=pooled[:, :, :-1])

    nb_alt = np.sum(
        np.apply_along_axis(
            lambda x: 1 if 1 in x else 0, axis=0, arr=pooled[:, :, :-1]
        )
    )

    if np.isin(pooled, -1).any():
        x = np.ones((1, 16, 1))
        y = np.asarray([-1, -1, 0])
        b = np.broadcast(x, y)
        decoded_gt = np.empty(b.shape)
        decoded_gt.flat = [u * v for (u, v) in b]

    else:
        encoded = np.dot(scores, d).reshape(1, 16, 1)
        b = np.broadcast_to(encoded, (1, 16, 3))
        if nb_alt == 0:
            decoded_gt = np.zeros_like(call)
        elif nb_alt == 2:
            decoder = lambda x: [1, -1, 0] if np.all(x == 2) else [0, 0, 0] # np.all() because of b.shape
            decoded_gt = np.apply_along_axis(decoder, axis=2, arr=b)
        else: # nb_alt > 2 and nb_alt < 8: # nb_alt = 2*n with n > 1
            decoder = lambda x: ([-1, -1, 0] if np.all(x == 2)
                                 else ([0, 0, 0] if (np.all(x == 1) or np.all(x == 0))
                                       else [1, 1, 0]))
            decoded_gt = np.apply_along_axis(decoder, axis=-1, arr=b)
    return decoded_gt


def _design_based(M):
    oo = pool.SNPsPool()
    d = oo.design_matrix()
    call = M['in'].reshape((1, 16, 3))

    pools_list = np.where(d == 1)[1].reshape((8,4))

    Z = decode_genotypes(d, call)

    return Z.reshape((4, 4, 3))


def test_pool_decode():
    for i, m in enumerate([M0, M1, M2, M3, M4, M5, M6, M7, M8, M9]):
    # for i, m in enumerate([M1]):
        z = _pool_decode(m)
        t = np.all(m['out'] == z)
        print(str(i) + ' --> ', t)
        if not t:
            print(z)
        print('')


def test_design_based():
    for i, m in enumerate([M0, M1, M2, M3, M4, M5, M6, M7, M8, M9]):
    #for i, m in enumerate([M9]):
        z = _design_based(m)
        t = np.all(m['out'] == z)
        print(str(i) + ' --> ', t)
        if not t:
            print(z)
        print('')

def idx_subsampling():
    samples = VCF('ALL.chr20.snps.gt.vcf.gz').samples
    subset = np.random.choice(samples, size=(16,), replace=False)
    p = np.argwhere(np.isin(np.asarray(samples), subset))

    t = np.all(np.asarray(samples)[p].sort() == subset.sort())
    print('Identical index sets? ', t)

    for v in VCF('ALL.chr20.snps.gt.vcf.gz'):
        print(v)
        pooled_samples = np.asarray(v.genotypes)
        break
    print('Initial: ', pooled_samples[p])

    oo = pool.SNPsPool()
    d = oo.design_matrix()
    call = M4['in'].reshape((1, 16, 3))
    dcd_gt = decode_genotypes(d, call, drop=False)
    print(dcd_gt)
    np.put_along_axis(pooled_samples, p, dcd_gt, axis=0)
    np.put_along_axis(pooled_samples, p, np.asarray([-1, -1, 0]), axis=0)
    print('Final: ', pooled_samples[p])


def cyvcf_threads():
    thr = VCF('ALL.chr20.snps.gt.vcf.gz', threads=500)
    thr.set_threads(500)
    frame = inspect.currentframe()
    print(inspect.getargvalues(frame))
    for i,t in enumerate(thr):
        print(t)
        if i == 5:
            break


def test_het_computing():
    for i, m in enumerate([M0, M1, M2, M3, M4, M5, M6, M7, M8, M9]):
        h = alltls.per_site_heteroz(gt_array=m['in'])
        print(str(i) + ' --> ', h)


def test_gtgl_converter():
    for i, m in enumerate([M0, M1, M2, M3, M4, M5, M6, M7, M8, M9]):
        print(str(i))
        g_out = np.apply_along_axis(alltls.map_gt_gl,
                                    -1,
                                    m['out'].reshape(16, 1, 3))
        print(g_out)


if __name__ == '__main__':
    os.chdir(prm.WD)
    # test_pool_decode()
    # test_design_based()
    # idx_subsampling()
    # cyvcf_threads()
    # test_het_computing()
    test_gtgl_converter()