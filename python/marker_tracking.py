from cyvcf2 import VCF
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import *
import itertools
from scripts.VCFPooling.poolSNPs import parameters as prm
from scripts.VCFPooling.poolSNPs import pool
from scripts.VCFPooling.poolSNPs.alleles import alleles_tools as alltls
from persotools.files import *

print('Load parameters'.ljust(80, '.'))
os.chdir(prm.WD)
plots_path = prm.PLOTS_PATH
plt.rcParams["figure.figsize"] = [10, 10]
plt.rcParams['figure.constrained_layout.use'] = True
f = 'ALL.chr20.snps.gt.vcf.gz'

print('Load raw variants'.ljust(80, '.'))
raw = VCF(prm.RAW['imp'])  # VCF iterator object
splits = pool.split_pools(VCF(f).samples, 16, seed=123)  # list of lists



"""
id          maf      maf_inter  err
rs6052163   0.442107 3          0.000000
rs1418029   0.443510 3          0.029167
rs2011489   0.448518 3          0.000000
rs2325843   0.459335 3          0.006250
rs11475696  0.486579 3          0.233333
rs202537    0.488982 3          0.000000
rs6054087   0.489383 3          0.118750
rs2281575   0.494992 3          0.000000
rs6042207   0.538462 3          0.000000
rs6139234   0.551082 3          0.000000
rs2326371   0.573518 3          0.006250
rs6417627   0.586939 3          0.002083
rs4350821   0.600962 3          0.008333
rs912111    0.603165 3          0.316667
"""


def count_alt(arr): # counts carriers if on raw samples
    nb_alt = np.sum(
        np.apply_along_axis(
            lambda x: 1 if 1 in x else 0, axis=0, arr=arr
        )
    )
    return nb_alt


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


def marker_stats(data, idn, verbose=False):
    m = None
    for n, v in enumerate(data(idn)):
        print(v)
        m = v
        if verbose:
            print('Number of known samples: ', v.num_called)
            print('Number of unknown samples: ', v.num_unknown)
            print('Number of heterozygotes: ', v.num_het)
            print('Number of alt homozygotes: ', v.num_hom_alt)
        if n == 0:
            break
    return m


def raw_gt(groups, data, var):
    sets = []
    res = []
    for gp in groups[0][0:prm.pools_imp]:
        sets.append(pool.SNPsPool().set_subset(gp))
    for p in sets:
        p.set_line_values(data.samples, var)
        pgt = p.get_call().reshape((1, p.size, 3))
        res.append(pgt)
        # print(pgt)
        # break
    return res


def pooled_gt(groups, data, var):
    sets = []
    res = []
    for gp in groups[0][0:prm.pools_imp]:
        sets.append(pool.SNPsPool().set_subset(gp))
    for p in sets:
        p.set_line_values(data.samples, var)
        pgt = p.pool_genotypes()
        res.append(pgt)
        # print(pgt)
        # break
    return res


def decoded_gt(groups, data, var_in):
    sets = []
    res = []
    for gp in groups[0][0:prm.pools_imp]:
        sets.append(pool.SNPsPool().set_subset(gp))
    for p in sets:
        p.set_line_values(data.samples, var_in)
        var_out = p.decode_genotypes(np.asarray(var_in.genotypes))
        idx = np.argwhere(np.isin(data.samples, p))
        pgt = np.asarray(var_out)[idx].reshape((1, p.size, 3))
        res.append(pgt)
        # print(pgt)
        # break
    return res


def imputed_gt(groups, data, idn):
    var = marker_stats(data, idn)
    sets = []
    res = []
    for gp in groups[0][0:prm.pools_imp]:
        sets.append(pool.SNPsPool().set_subset(gp))
    for p in sets:
        p.set_line_values(data.samples, var)
        idx = np.argwhere(np.isin(data.samples, p.flatten()))
        pgt = np.asarray(var.genotypes)[idx].reshape((1, len(idx), 3))
        res.append(pgt)
    return res


def frame_comb_count():
    cb = itertools.combinations_with_replacement(range(1, 5), 2)
    comb = [(0, 0)] + [c for c in cb]
    df = pd.DataFrame.from_records(data=comb, columns=['N1', 'N2'])
    df['count'] = pd.Series(np.zeros((len(comb),), dtype=int))
    return df


def barplot_comb(df, lab, ttl):
    ax = df.plot.bar(title=ttl)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, lab,
              bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.show()


def boxplot_err_comb(df, cols=None):
    df.boxplot(by=['nbRow', 'nbCol'], column=cols)
    plt.legend
    plt.show()


if __name__ == '__main__':
    decoding_power()
    targets = {'pooled': {'id': ['rs1418029', 'rs2011489', 'rs2325843', 'rs11475696',
                                 'rs138183721', 'rs138276627', 'rs139656140'],
                          'pos': ['20:2060151', '20:888523', '20:2348882', '20:5749550',
                                  '20:3108058', '20:686062', '20:1169354'],
                          'maf': [0.443510, 0.448518, 0.459335, 0.486579,
                                  0.000200, 0.000200, 0.000200],
                          'err': [0.029167, 0.0, 0.006250, 0.233333,
                                  0.0, 0.152083, 0.002083]},
               'missing': ['rs73911387', '20:50359566']}  # rs6031147 for missing

    err = pd.read_csv('pooled.sorted.chunk{}.csv'.format(prm.CHK_SZ),
                      sep='\t',
                      header=[0, 1],
                      index_col=[0, 1, 2])
    # Sort by ascending MAF, descending error
    eps = err.mean(axis=1).rename('err').to_frame()
    slc = eps.sort_values(by='err', ascending=False).sort_index(axis=0,
                                                                level='maf',
                                                                ascending=True)
    print(slc)

    d1, d2 = list(), list()
    for it in range(len(targets['pooled']['pos'])):
        mrk = marker_stats(raw, targets['pooled']['pos'][it], verbose=True)
        g0 = raw_gt(splits, raw, mrk)
        g1 = pooled_gt(splits, raw, mrk)
        g2 = decoded_gt(splits, raw, mrk)
        g3 = imputed_gt(splits, VCF(prm.POOLED['corr'] + '.vcf.gz'), targets['pooled']['pos'][it])

        df0 = frame_comb_count()
        df1 = pd.DataFrame(data=np.empty((prm.pools_imp, 4),
                                         dtype=float),
                           columns=['nbRow', 'nbCol', 'missing_alleles', 'imp_err'])
        for i, g in enumerate(itertools.zip_longest(g0, g1, g2, g3)):
            df1.loc[i, 'nbRow'] = count_alt(g[1][:, :4, :-1])
            df1.loc[i, 'nbCol'] = count_alt(g[1][:, 4:, :-1])
            diff03 = alltls.para_array(g[0][:, :, :-1], g[3][:, :, :-1])
            df1.loc[i, 'imp_err'] = alltls.per_axis_error(diff03, 0)
            df1.loc[i, 'missing_alleles'] = np.sum(alltls.count_missing_alleles(gt_array=g[2]))
        df1 = df1.astype({'nbRow': int, 'nbCol': int}, copy=True)
        d1.append(df1)
        # boxplot_err_comb(df1, cols=['missing_alleles'])
        # boxplot_err_comb(df1, cols=['imp_err'])

        df2 = df1.groupby(['nbRow', 'nbCol'])['missing_alleles'].count()#.rename(columns={'missing_alleles': 'count'})
        d2.append(df2)

        for i, j in df2.index:
            tot_ij = df2.xs((i, j))
            if i != j:
                try:
                    tot_ij += df2.xs((j, i))
                except:
                    pass
            df0.loc[(df0['N1'] == i) & (df0['N2'] == j), 'count'] = tot_ij
        df0.reset_index(drop=True, inplace=True)
        df0.set_index(['N1', 'N2'], inplace=True)
        if it == 0:
            d0 = df0
        else:
            d0 = d0.join(df0, rsuffix='_' + targets['pooled']['id'][it])

    print(d0)
    barplot_comb(d0,
                 ttl='Markers: (id, maf, imputation error)',
                 lab=list(zip(targets['pooled']['id'],
                                  targets['pooled']['maf'],
                                  targets['pooled']['err'])))

"""
Select a marker with low MAF and high error in imputation:
1) Analyze patterns of pools (sort by nbRow, nbCol with ALT allele)
2) Compute decoding power based on designed matrix
2) Analyze imputation error in each pool and plot it on x-axis=nbRow, y-axis=nbCol
"""


