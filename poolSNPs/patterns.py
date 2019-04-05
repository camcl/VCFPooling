from cyvcf2 import VCF
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import *
import itertools
from scripts.poolSNPs import parameters as prm
from scripts.poolSNPs import pool
from scripts.poolSNPs.alleles import alleles_tools as alltls
from scripts.poolSNPs.alleles import alleles_plots as allplt
from scripts.poolSNPs import results as res
from persotools.debugging import *
from persotools.files import *


def combine_gl(nbSlots, nbIdv, gtpos=['RA', 'AA']):
    prd = itertools.combinations_with_replacement(gtpos, r=nbIdv)
    patterns = list()
    for p in prd:
        tpl = p + tuple(['RR'] * (nbSlots - nbIdv))
        patterns.append(tpl)
    pot = np.array(patterns).flatten()
    # adjusted GLs
    wgl = np.array([np.count_nonzero(np.where(pot == let, 1, 0)) for let in ['RR', 'RA', 'AA']]) / len(pot)
    return wgl.squeeze()


def generate_patterns(nbRow, nbCol, nbIdv):
    n = max((nbRow, nbCol))
    N = nbRow * nbCol
    patset = list()
    comb = itertools.combinations(range(N), r=nbIdv)
    for cb in comb:
        p = np.isin(np.arange(0, N), cb).reshape((nbRow, nbCol))
        patset.append(p.astype(int))
    return patset


def count_pos(arr): # counts carriers if on raw samples
    nbPos = lambda x: 1 if 1 in x else 0
    pattern = np.where(np.isin(arr, -1), 1, 0)
    if len(pattern) == 0:
        return 0
    else:
        rowPos = np.apply_along_axis(
                nbPos, axis=0, arr=pattern
            )
        colPos = np.apply_along_axis(
            nbPos, axis=1, arr=pattern
        )
        return np.sum(rowPos)//2, np.sum(colPos)//2


def count_pattern(arr): # counts carriers if on raw samples
    nbPos = lambda x: 1 if 1 in x else 0
    rowPos = np.apply_along_axis(
            nbPos, axis=0, arr=arr
        )
    colPos = np.apply_along_axis(
        nbPos, axis=1, arr=arr
    )
    return np.sum(rowPos), np.sum(colPos)


def constraint_pattern(nbRow, nbCol, nbIdv):
    patset = generate_patterns(nbRow, nbCol, nbIdv)
    res = 0
    for p in patset:
        # print(p)
        r, c = count_pattern(p)
        if r * c < np.size(p):
            continue
            # print('Reject')
        else:
            #print('Accept')
            res += 1
        # print('')
    return res


def screen_patterns(nbRow, nbCol):
    n = max((nbRow, nbCol))
    N = nbRow * nbCol
    for i in range(n, N+1, 1):
        tot = constraint_pattern(nbRow, nbCol, i)
        print(i)
        print('total --> ', tot)


def weight_gl(nbRow, nbCol):
    if nbRow < 2 or nbCol < 2:
        w = np.array([0, 0.5, 0.5])
    else:
        n = max((nbRow, nbCol))
        N = nbRow * nbCol
        gl = np.zeros(3)
        den = 0
        for i in range(n, N + 1, 1):
            tot = constraint_pattern(nbRow, nbCol, i)
            ratios = combine_gl(N, i)
            gl = np.sum(np.vstack([gl, tot * ratios]), axis=0)
            den += tot
        w = gl / den
    return w


def enumerate_comb():
    cb = itertools.combinations_with_replacement(range(1, 4+1), 2)
    weighted = {('0', '0'): np.array([1/3, 1/3, 1/3])}
    for c in cb:
        weighted[(str(c[0]), str(c[1]))] = weight_gl(c[0], c[1])
    return weighted


def decoding_power(mtx=pool.SNPsPool().design_matrix()):
    """
    Weight w = number of pools in which a specimen is participating.
    Intersection lbd = dot-product of 2 column vectors.
    Reveals the number of intersecting pools that hold 2 particular specimens.
    d = (w_min -1)/lbd_max
    :param mtx: design matrix. Describes the pooling pattern used.
    :return: d = decoding power
    """
    w = np.sum(mtx, axis=1)
    cicj = itertools.combinations(range(mtx.shape[1]), 2)
    lbd = np.array([np.transpose(mtx[:,i]).dot(mtx[:,j]) for i,j in cicj])
    d = (np.min(w) - 1)/np.max(lbd)
    # print(w)
    # print(lbd)
    print('decoding power = ', int(d))
    return d


def replace_gl(variant, d):
    # reshaping as square-looking pools
    g = np.array(variant.genotypes)
    g_in = g.reshape(
        (len(g) // prm.pools_size,
         prm.pools_size,
         3)
    )
    g_out = list()

    for p in g_in:
        nbRow, nbCol = count_pos(p.reshape(4, 4, 3)[:, :, :-1])
        r = min(nbRow, nbCol)
        c = max(nbRow, nbCol)
        gl = d[(str(r), str(c))]
        t = lambda x: alltls.map_gt_gl(x, unknown=gl)
        q = np.apply_along_axis(t, axis=-1, arr=p)
        g_out.append(q)
    g_out = np.array(g_out).reshape(g.shape)
    return g_out


def adaptative_likelihood_converter(f_in, f_out):
    d = enumerate_comb()
    replace = lambda v: replace_gl(v, d)
    alltls.file_likelihood_converter(f_in,
                                     f_out,
                                     func=replace)


# w22 = weight_gl(2, 2)
# print('2x2\r')
# print(w22)
# print('\r\n')
# w23 = weight_gl(2, 3)
# print('2x3\r')
# print(w23)
# print('\r\n')
# w33 = weight_gl(3, 3)
# print('3x3\r')
# print(w33)
# print('\r\n')
# w34 = weight_gl(3, 4)
# print('3x4\r')
# print(w34)
# print('\r\n')
# w44 = weight_gl(4, 4)
# print('4x4\r')
# print(w44)

if __name__ == '__main__':
    os.chdir(prm.WD)
    d = enumerate_comb()
    print(d)
    replace = lambda v: replace_gl(v, d)
    alltls.file_likelihood_converter('IMP.chr20.pooled.snps.gt.chunk1000.vcf.gz',
                                     'IMP.chr20.pooled.snps.adaptative_gl.chunk1000.vcf',
                                     func=replace)
