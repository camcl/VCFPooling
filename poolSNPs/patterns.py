import numpy as np
import pandas as pd
from scipy.stats import *
import itertools
import multiprocessing as mp
from collections import Counter

from scripts.poolSNPs import parameters as prm
from scripts.poolSNPs import pool
from scripts.poolSNPs.alleles import alleles_tools as alltls
from scripts.poolSNPs.alleles import alleles_plots as allplt
from scripts.poolSNPs import results as res
from persotools.debugging import *
from persotools.files import *

import inspect
import time


def report(xk):
    frame = inspect.currentframe().f_back
    #print(frame)
    locals = frame.f_locals
    #print(locals)
    return locals


def combine_gt(nbIdv, gtpos=['RA', 'AA']):
    prd = itertools.combinations_with_replacement(gtpos, r=nbIdv)
    return prd


def permute_gt(nbSlots, p):
    tpl = p + tuple(['RR'] * (nbSlots - len(p)))
    t = itertools.permutations(tpl)
    iter_set = set(s for s in t).__iter__()
    return iter_set


def parallel_permut(nbSlots, prd):
    ptask = mp.Pool(processes=4)
    args = tuple(zip(itertools.repeat(nbSlots),
                     [p for p in prd]))
    results = ptask.starmap(permute_gt, args)
    ptask.close()
    ptask.join()
    patterns = iter(r for r in results)
    return patterns


def generate_patterns(nbRow, nbCol, nbIdv):
    combi = combine_gt(nbIdv)
    patterns = parallel_permut(nbRow*nbCol, combi)
    patset = list()
    t = np.vectorize(lambda x: 0 if x == 'RR' else 1)
    for p in patterns:
        for a in p:
            q = np.apply_along_axis(t,
                                    axis=0,
                                    arr=a)
            patset.append((a, q))
    return patset.__iter__()


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
    # filters patterns
    patin = generate_patterns(nbRow, nbCol, nbIdv)
    patout = list()
    for p in patin:
        r, c = count_pattern(p[1].reshape((nbRow, nbCol)))
        if r * c < np.size(p[1]):
            continue
            # print('Reject')
        else:
            #print('Accept')
            patout.append(p[0])
    return itertools.chain.from_iterable(patout)


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
        pot = itertools.starmap(constraint_pattern,
                                zip(itertools.repeat(nbRow, times=N-n+1),
                                    itertools.repeat(nbCol, times=N-n+1),
                                    range(n, N+1))
                                )
        cnt = Counter(itertools.chain.from_iterable(pot))
        ratios = np.array([cnt['RR'], cnt['RA'], cnt['AA']])
        w = ratios / sum(([cnt['RR'], cnt['RA'], cnt['AA']]))
        print('Processed: ({}, {}) {}'.format(nbRow, nbCol, w))
    return w


def enumerate_comb():
    cb = itertools.combinations_with_replacement(range(1, 4+1), 2)
    weighted = {('0', '0'): np.array([1/3, 1/3, 1/3])}
    for c in cb:
        weighted[(str(c[0]), str(c[1]))] = weight_gl(c[0], c[1])
        # save in a file!
    df = pd.DataFrame.from_dict(data=list([k[0], k[1], i] for k,i in weighted.items()),
                                orient='index',
                                columns=['n', 'm', 'gl'])
    df.to_csv(os.path.join(prm.WD, 'adaptative_gl.csv'))


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


def adaptative_likelihood_converter(f_in, f_out, write=False):
    if write:
        # write gl and combinations
        enumerate_comb()
    df = pd.read_csv(os.path.join(prm.WD, 'adaptative_gl.csv'))
    d = df.to_dict()
    replace = lambda v: replace_gl(v, d)
    alltls.file_likelihood_converter(f_in,
                                     f_out,
                                     func=replace)


# if __name__ == '__main__':
os.chdir(prm.WD)

start = time.time()
enumerate_comb()
stop = time.time()
print('Elapsed time: ', stop-start)
