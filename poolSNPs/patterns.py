import pandas as pd
from scipy.stats import *
import itertools
import multiprocessing as mp
from collections import Counter
import functools
import operator
import math

from cyvcf2 import VCF
from scripts.VCFPooling.poolSNPs import parameters as prm
from scripts.VCFPooling.poolSNPs import pool
from scripts.VCFPooling.poolSNPs.alleles import alleles_tools as alltls
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
    print('\nnbIdv:', nbIdv)
    combi = combine_gt(nbIdv)
    n = max((nbRow, nbCol))
    # squeeze patterns to test using symmetric behavior of RA and AA
    if nbIdv - n >= 0 and (nbIdv - n)%2 == 0:
        try:
            combi = combi[:(nbIdv-n+1)//2]
        except TypeError:
            pass
    if nbIdv - n >= 0 and (nbIdv - n)%2 != 0:
        try:
            combi = combi[:1 + (nbIdv-n+1)//2]
        except TypeError:
            pass
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


def count_pos(arr, pos=1): # counts carriers if on raw samples
    nbPos = lambda x: 1 if 1 in x else 0
    pattern = np.apply_along_axis(lambda x: 1 if pos in x else 0,
                                  -1,
                                  arr)
    rowPos = np.apply_along_axis(
        nbPos, axis=0, arr=pattern
    )
    colPos = np.apply_along_axis(
        nbPos, axis=1, arr=pattern
    )
    return np.sum(rowPos), np.sum(colPos)


def count_pattern(arr): # counts carriers if on translated patterns
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
    cnt = {'AA': 0, 'RA': 0, 'RR': 0}
    cnt = {**cnt, **dict(Counter(itertools.chain.from_iterable(patout)))}
    print(cnt)
    return cnt


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
                                zip(itertools.repeat(nbRow, times=N-n),
                                    itertools.repeat(nbCol, times=N-n),
                                    range(n, N))
                                )
        cnt = dict(functools.reduce(operator.add, map(Counter, pot)))
        ratios = np.array([cnt['RR']*2, cnt['RA'] + cnt['AA'], cnt['AA'] + cnt['RA']])
        w = ratios / sum(ratios)
        print('Processed: ({}, {}) {}'.format(nbRow, nbCol, w))
    return w


def enumerate_comb():
    cb = itertools.combinations_with_replacement(range(1, 4), 2)
    weighted = {('0', '0'): np.array([1/3, 1/3, 1/3])}
    for c in cb:
        weighted[(str(c[0]), str(c[1]))] = weight_gl(c[0], c[1])
        # save in a file!
    df = pd.DataFrame.from_dict(data=list([k[0], k[1], i] for k, i in weighted.items()),
                                orient='index',
                                columns=['n', 'm', 'rr', 'ra', 'aa'])
    df.to_csv(os.path.join(prm.WD, 'adaptative_gl.csv'), header=False, index=False)


def decoding_power():
    """
    Weight w = number of pools in which a specimen is participating.
    Intersection lbd = dot-product of 2 column vectors.
    Reveals the number of intersecting pools that hold 2 particular specimens.
    d = (w_min -1)/lbd_max
    :param mtx: design matrix. Describes the pooling pattern used.
    :return: d = decoding power
    """
    mtx = pool.SNPsPool().design_matrix()
    w = np.sum(mtx, axis=1)
    cicj = itertools.combinations(range(mtx.shape[1]), 2)
    lbd = np.array([np.transpose(mtx[:,i]).dot(mtx[:,j]) for i,j in cicj])
    d = (np.min(w) - 1)/np.max(lbd)
    print('decoding power = ', int(d))
    return d


def replace_gl(variant, d, log=True):
    # reshaping as square-looking pools
    g = np.array(variant.genotypes)
    g_in = g.reshape(
        (len(g) // prm.pools_size,
         prm.pools_size,
         3)
    )
    g_out = list()

    convert = lambda x: alltls.map_gt_gl(x, unknown=gl)
    logzero = np.vectorize(lambda x: -5.0 if x <= pow(10, -5) else math.log10(x))

    for p in g_in:
        nbRow, nbCol = count_pos(p.reshape((4, 4, 3))[:, :, :-1], pos=-1)
        r = min(nbRow, nbCol)
        c = max(nbRow, nbCol)
        if r * c == 0:
            q = p
        else:
            gl = d[(r, c)]
            q = np.apply_along_axis(convert, axis=-1, arr=p)
        g_out.append(q)
    g_out = np.array(g_out).reshape(g.shape).astype(float)

    if log:
        g_out = logzero(g_out)
    return g_out


def adaptative_likelihood_converter(f_in, f_out, write=False):
    if write:
        # write gl and combinations
        enumerate_comb()
    df = pd.read_csv(os.path.join(prm.WD, 'adaptative_gl.csv'),
                     header=None,
                     names=['n', 'm', 'rr', 'ra', 'aa']
                     )
    df2dict = dict(((int(n), int(m)), [rr, ra, aa]) for n, m, rr, ra, aa in df.itertuples(index=False, name=None))
    replace = lambda v: replace_gl(v, df2dict)
    print('f_in --> ', f_in)
    print('f_out --> ', f_out)
    alltls.file_likelihood_converter(f_in,
                                     f_out,
                                     func=replace)


def permute_4x4(nbIdv, p_list):
    # number of combinations of 2 items among 12 ones
    nbplac = list(itertools.combinations(range(12), nbIdv))
    ini = {'AA': 0, 'RA': 0, 'RR': (12-nbIdv) * len(p_list)}
    cnt = dict(Counter(itertools.chain.from_iterable(p_list)))
    cnt = {**ini, **cnt}
    cnt = dict((k, v * len(nbplac)) for k,v in cnt.items())
    # 16 possible patterns for the min pattern...
    expand = np.array([cnt['RR'], cnt['RA'], cnt['AA']]) * 16
    # ... adding 40 RA and 40 AA...
    expand = np.add(expand, [0, 40, 40])
    # ... giving the result for only 1 side of the symmetry
    expand = expand * 2
    cnt = dict(zip(['RR', 'RA', 'AA'], expand))
    print(cnt)
    return cnt


def generate_patterns_4x4(nbIdv, n=4):
    print('\nnbIdv:', nbIdv)
    combi = combine_gt(nbIdv)
    patterns = permute_4x4(nbIdv, list(combi))
    return patterns


def case_4x4(n=4, N=12, write=False):
    # minpat = np.identity(n).reshape((n**2,))
    pot = map(generate_patterns_4x4, range(N+1))
    cnt = dict(functools.reduce(operator.add, map(Counter, pot)))
    ratios = np.array([cnt['RR'], cnt['RA'], cnt['AA']])
    w = ratios / sum(ratios)
    print([n, n] + w.tolist())
    if write:
        df = pd.read_csv(os.path.join(prm.WD, 'adaptative_gl.csv'),
                         header=None,
                         names=['n', 'm', 'rr', 'ra', 'aa'])
        df = df.append(pd.Series([n, n] + w.tolist(), index=['n', 'm', 'rr', 'ra', 'aa']),
                       ignore_index=True
                       )
        df2dict = dict(((int(n), int(m)), [rr, ra, aa]) for n, m, rr, ra, aa in df.itertuples(index=False, name=None))
        print(df2dict)
        df.to_csv(os.path.join(prm.WD, 'adaptative_gl.csv'), header=False, index=False)
    return w


if __name__ == '__main__':
    os.chdir(prm.WD)
    start = time.time()
    # enumerate_comb()
    # w44 = case_4x4(write=False)
    # stop = time.time()
    # print('\r\n', w44)
    # print('Elapsed time: ', stop-start)
    df = pd.read_csv(os.path.join(prm.WD, 'adaptative_gl.csv'),
                     header=None,
                     names=['n', 'm', 'rr', 'ra', 'aa'])
    df2dict = dict(((int(n), int(m)), [rr, ra, aa]) for n, m, rr, ra, aa in df.itertuples(index=False, name=None))

    data = VCF(os.path.join(prm.WD, 'gt', 'IMP.chr20.pooled.snps.gt.chunk1000.vcf.gz'))
    for i, var in enumerate(data('20:51044694-62887012')):
        v = var
        if i == 1:
            break
        g_in = np.array(v.genotypes)
        g_out = replace_gl(v, df2dict)
        z = list(zip(g_in, g_out))
        for pairz in z:
            print(pairz[0], pairz[1])
        print('\n')

