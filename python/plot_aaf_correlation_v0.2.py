import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from cyvcf2 import VCF

from scripts.poolSNPs.alleles import alleles_tools as alltls
from scripts.poolSNPs.alleles import alleles_plots as allplt
from scripts.poolSNPs import parameters as prm
from persotools.files import *

"""
Plot correlation between AAF after imputation (computed with Phaser)
and AAF in the original dataset (theoretical).
Compare effects of imputation:
- GT vs GL
- varying GL values for filling in missing genotypes
"""


def get_df_evol(wd: str, k_set: list, path_all: str, idt='id') -> pd.DataFrame:
    """
    Compute alternate alleles frequencies after before/after imputation
    :param wd:
    :param keys:
    :param path_all:
    :param idt:
    :return:
    """
    df_aaf = alltls.get_aaf(path_all, idt=idt)

    list_aaf = alltls.compute_aaf_evol(wd, k_set, idt=idt)
    df_aaf = df_aaf.join(list_aaf)
    df_aaf.sort_values(by='af_info', axis=0, inplace=True)

    return df_aaf


def plot_aaf_correlation(df_aaf, col_set, typ='scatter'):
    plt.rcParams["figure.figsize"] = [12, 6]
    plt.rcParams["figure.autolayout"] = True

    lin_aaf, ax_lin = plt.subplots()
    for a, k_set in enumerate(col_set):
        colors = tuple(np.random.rand(3))  # RGB --> 3
        if typ == 'scatter':
            df_aaf.plot.scatter(x='af_info',
                                y='preimp_' + k_set,
                                ax=ax_lin,
                                label='preimp_' + k_set,
                                marker='o',
                                color=np.concatenate((colors, [0.2])))  # RGBA, transparency
            df_aaf.plot.scatter(x='af_info',
                                y='postimp_' + k_set,
                                ax=ax_lin,
                                label='postimp_' + k_set,
                                marker='o',
                                color=colors)
        plt.xlabel('Theoretical AAF')
        plt.ylabel('Dataset AAF')
        plt.plot(range(2), linestyle='-', color='k')
        plt.legend()

    plt.title('AAF evolution through processing of data sets', loc='center')
    plt.suptitle("")
    plt.savefig('aaf_evol.{}.gtgl.pooled.chunk{}.png'.format(typ, df_aaf.shape[0]),
                orientation='landscape',
                dpi='figure')


if __name__ == '__main__':
    # Configure working directory
    print('Configure working directory'.ljust(80, '.'))
    if prm.GTGL == 'GT':
        cd = prm.PATH_GT_FILES
    if prm.GTGL == 'GL':
        if prm.unknown_gl != 'adaptative':
            cd = os.path.join(prm.WD, 'gl', 'gl_' + '_'.join(np.around(prm.unknown_gl, decimals=2).astype(str)))
        else:
            cd = os.path.join(prm.WD, 'gl', 'gl_adaptive', 'phaser')
    os.chdir(cd)

    print('Load parameters'.ljust(80, '.'))
    sorting = True  # sort data sets by AAF and population values
    params = [('gt', prm.CHK_SZ, 'id'), ('gl', prm.CHK_SZ, 'id')]
    devol = []

    # Load AAFs
    allaafs = alltls.get_aaf(prm.PATH_GT_FILES + '/ALL.chr20.pooled.snps.gt.chunk{}.vcf.gz'.format(prm.CHK_SZ),
                             idt='id')
    aafs = allaafs.loc[:, 'af_info'].to_frame()
    impaafs = alltls.get_aaf(prm.PATH_GT_FILES + '/IMP.chr20.pooled.snps.gt.chunk{}.vcf.gz'.format(prm.CHK_SZ),
                             idt='id')
    compaafs = allaafs.loc[:, ['af_info', 'aaf']].join(impaafs.loc[:, ['af_info', 'aaf']], how='inner', rsuffix='_imp')

    for p_set in params:
        print('\nSet params:', p_set)
        gtgl, sz, idt = p_set
        ALL = prm.PATH_GT_FILES + '/ALL.chr20.snps.gt.chunk{}.vcf.gz'.format(sz)
        B1 = prm.PATH_GT_FILES + '/IMP.chr20.pooled.snps.gt.chunk{}.vcf.gz'.format(sz)
        if gtgl == 'gl':
            POOL = os.path.join(cd,
                                'IMP.chr20.pooled.imputed.gt.chunk{}.vcf.gz'.format(sz))
        if gtgl == 'gt':
            POOL = os.path.join(prm.PATH_GT_FILES,
                                'all_snps_all_samples',
                                'IMP.chr20.pooled.imputed.gt.chunk{}.vcf.gz'.format(sz))

        df_chk = get_df_evol(os.path.dirname(POOL),
                             'pooled',
                             ALL,
                             idt=idt).loc[:, ['preimp_pooled', 'postimp_pooled']]
        df_chk.rename(columns={'preimp_pooled': 'preimp_' + str(gtgl),
                               'postimp_pooled': 'postimp_' + str(gtgl)},
                      inplace=True)

        devol.append(df_chk)

    df_plot = aafs.join(devol, how='inner')
    df_plot.sort_values(by='af_info', inplace=True)
    # print('dfplot\n', df_plot)
    plot_aaf_correlation(df_plot, col_set=list(zip(*params))[0], typ='scatter')
    df_err = pd.read_csv('pooled' + '{}.chunk{}.csv'.format('.sorted' if sorting else '', prm.CHK_SZ),
                         sep='\t',
                         encoding='utf-8',
                         index_col=0,
                         usecols=[0] + list(range(3, prm.NB_IMP + 3)),
                         skiprows=[0, 1]
                         )
    allplt.plot_aaf_twist('pooled',
                          df_err,
                          os.path.join(cd, 'IMP.chr20.pooled.imputed.gt.chunk{}.vcf.gz'.format(prm.CHK_SZ)))

