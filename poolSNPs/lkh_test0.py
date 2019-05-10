import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from cyvcf2 import VCF

from scripts.poolSNPs.alleles import alleles_tools as alltls
import scripts.poolSNPs.results as res
from scripts.poolSNPs import parameters as prm
from persotools.files import *
from persotools.debugging import *

"""
Convert GT to GL before imputing.
Objective: get same plot when not modifying ?-GL:
"""


def compute_maf_evol(set, fpre, fpost, gtgl='gt', chk_sz=None, idt='id'):
    """

    :param set: str, short name of the data set pooled/missing...
    :param df_maf: maf from the original vcf-chunked file, with markers ID as index
    :return:
    """
    if chk_sz is None:
        chk_sz = prm.CHK_SZ
    print('\r\nSet --> ', set)
    steps = {'preimp': fpre,
             'postimp': fpost}
    temp = []
    for d, vcf in steps.items():
        df = pd.DataFrame.from_dict(alltls.compute_maf(vcf, idt=idt),
                                    orient='index',
                                    columns=[d + '_' + set])
        temp.append(df)

    return temp


def get_df_evol(err_dic, fpre, fpost, path_all, gen, chk, idt='id'):
    bin_maf = prm.BIN_MAF

    df_maf = alltls.get_maf(path_all, id=idt)

    for k_set in err_dic.keys():
        list_maf = compute_maf_evol(k_set, fpre, fpost, gtgl=gen, chk_sz=chk, idt=idt)
        df_maf = df_maf.join(list_maf)
    df_maf.sort_values(by='maf', axis=0, inplace=True)

    if bin_maf:
        convert = np.vectorize(lambda x: alltls.convert_maf(x))
        inter = prm.INTER
        np_bin_maf = np.digitize(convert(df_maf.values),
                                 bins=inter)
        np_bin_maf = np.subtract(np_bin_maf/len(inter), 0.00)
        df_maf = pd.DataFrame(data=np_bin_maf, index=df_maf.index, columns=df_maf.columns)

    return df_maf


def plot_test(df_maf, col_set, typ='scatter'):
    plt.rcParams["figure.figsize"] = [12, 6]
    plt.rcParams["figure.autolayout"] = True

    lin_maf, ax_lin = plt.subplots()
    for a, k_set in enumerate(col_set):
        colors = tuple(np.random.rand(3)) # RGBA --> 3
        if typ == 'scatter':
            df_maf.plot.scatter(x='maf',
                                y='preimp_' + k_set,
                                ax=ax_lin,
                                label='preimp_' + k_set,
                                marker='v',
                                color=colors)
            df_maf.plot.scatter(x='maf',
                                y='postimp_' + k_set,
                                ax=ax_lin,
                                label='postimp_' + k_set,
                                marker='o',
                                color=colors)
        if typ == 'line':
            df_maf.plot(x='maf',
                        y='preimp_' + k_set,
                        ax=ax_lin,
                        linestyle='--',
                        color=colors)
            df_maf.plot(x='maf',
                        y='postimp_' + k_set,
                        ax=ax_lin,
                        linestyle='-',
                        color=colors)
        plt.xlabel('Theoretical MAF')
        plt.ylabel('Dataset MAF')
        plt.plot(range(2), linestyle='-', color='k')
        delta_post = df_maf['postimp_' + k_set].sub(df_maf['maf'])
        err_post = alltls.rmse_df(delta_post.to_frame(), kind='mse')
        plt.text(0.65, 0.25-0.1*a, 'MSE postimp_' + k_set + ' = ' + str(err_post))
        plt.legend()

    plt.title('MAF evolution through processing of data sets', loc='center')
    plt.suptitle("")
    plt.savefig('maf_evol.{}.gtgl.pooled.chunk{}.png'.format(typ, df_maf.shape[0]),
                orientation='landscape',
                dpi='figure')


if __name__ == '__main__':
    ### DATA
    if prm.GTGL == 'GT':
        cd = os.path.join(prm.WD, 'gt/stratified')
    if prm.GTGL == 'GL':
        if prm.unknown_gl != 'adaptative':
            cd = os.path.join(prm.WD, 'gl', 'gl_' + '_'.join(np.around(prm.unknown_gl, decimals=2).astype(str)))
        else:
            # cd = os.path.join(prm.WD, 'gl')
            cd = os.path.join(prm.WD, 'gl', 'gl_adaptative') # , 'kristiina')
    print('Load parameters'.ljust(80, '.'))
    sorting = True # sort data sets by MAF and population values
    os.chdir(cd)
    params = [('gt', prm.CHK_SZ, 'chrom:pos')]#, ('gl', prm.CHK_SZ, 'chrom:pos')]
    devol = []

    # Load MAFs
    allmafs = alltls.get_maf(prm.WD + '/gt/stratified/ALL.chr20.pooled.snps.gt.chunk{}.vcf.gz'.format(prm.CHK_SZ),
                             id='chrom:pos')
    mafs = allmafs.loc[:, 'maf'].to_frame()
    aafs = allmafs.loc[:, 'aaf'].to_frame()
    impmafs = alltls.get_maf(prm.WD + '/gt/stratified/IMP.chr20.pooled.snps.gt.chunk{}.vcf.gz'.format(prm.CHK_SZ),
                             id='chrom:pos')
    compmafs = allmafs.loc[:, ['maf', 'aaf']].join(impmafs.loc[:, ['maf', 'aaf']], how='inner', rsuffix='_imp')

    for p_set in params:
        print('\nSet params:', p_set)
        gtgl, sz, idt = p_set
        ALL = prm.WD + '/gt/stratified/ALL.chr20.snps.gt.chunk{}.vcf.gz'.format(sz)
        B1 = prm.WD + '/gt/stratified/IMP.chr20.pooled.snps.gt.chunk{}.vcf.gz'.format(sz)
        if gtgl == 'gl':
            POOL = 'IMP.chr20.pooled.imputed.gt.chunk{}.vcf.gz'.format(sz)
        if gtgl == 'gt':
            POOL = prm.WD + '/gt/stratified/IMP.chr20.pooled.imputed.gt.chunk{}.vcf.gz'.format(sz)

        print('Raw samples and index labels'.ljust(80, '.'))
        # Create line-index and column-index (multi-indexing)
        maf_idx, pop_idx = res.make_index(B1, src=ALL, id=idt)

        ### PROCESSING
        print('Build datasets'.ljust(80, '.'))
        raw0, raw1 = alltls.vcf2dframe(VCF(B1), prm.NB_IMP, id=idt)
        pool0, pool1 = alltls.vcf2dframe(VCF(POOL), prm.NB_IMP, id=idt)
        # raw0 = raw0.join(pool0.index.to_frame(), how='inner').iloc[:, :-1]
        # raw1 = raw1.join(pool1.index.to_frame(), how='inner').iloc[:, :-1]

        raw0 = raw0.join(mafs, how='inner').drop('maf', axis=1)
        raw1 = raw1.join(mafs, how='inner').drop('maf', axis=1)
        pool0 = pool0.join(mafs, how='inner').drop('maf', axis=1)
        pool1 = pool1.join(mafs, how='inner').drop('maf', axis=1)

        samples = VCF(B1).samples
        variants = raw0.index.tolist()

        # NumPy juxtapos for error computation
        raw = np.stack([raw0.values, raw1.values], axis=-1)
        pool = np.stack([pool0.values, pool1.values], axis=-1)

        set_errors = alltls.compute_imp_err(['pooled'],
                                            [pool],
                                            raw,
                                            maf_idx,
                                            pop_idx,
                                            sorting)
        df_chk = get_df_evol(set_errors,
                             B1,
                             POOL,
                             ALL,
                             gtgl,
                             sz,
                             idt=idt).loc[:, ['preimp_pooled', 'postimp_pooled']]
        df_chk.rename(columns={'preimp_pooled': 'preimp_' + str(gtgl),
                               'postimp_pooled': 'postimp_' + str(gtgl)},
                      inplace=True)

        devol.append(df_chk)

    df_plot = mafs.join(devol, how='inner')
    print(df_plot)
    df_plot.sort_values(by='maf', inplace=True)
    print('dfplot\n', df_plot)
    plot_test(df_plot, col_set=list(zip(*params))[0], typ='line')
    plot_test(df_plot, col_set=list(zip(*params))[0], typ='scatter')

