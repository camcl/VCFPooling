import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from itertools import starmap
from functools import partial
from cyvcf2 import VCF

from scripts.VCFPooling.python.archived.alleles import alleles_tools as alltls
from scripts.VCFPooling.poolSNPs import dataframe as vcfdf
from scripts.VCFPooling.poolSNPs import parameters as prm
from persotools.files import *

"""
Plot correlation between AAF after imputation (computed)
and AAF in the original dataset (theoretical).
Compare effects of imputation:
- GT vs GL
- completely remove 121 markers separately and concatenate them together with the remaining for getting 10000
 v.s impute 10000 markers in 1 batch
"""


def compute_aaf_evol(set, fpre, fpost, gtgl='gt', chk_sz=None, idt='id'):
    """

    :param set: str, short name of the data set pooled/missing...
    :param df_aaf: aaf from the original vcf-chunked file, with markers ID as index
    :return:
    """
    if chk_sz is None:
        chk_sz = prm.CHK_SZ
    print('\r\nSet --> ', set)
    steps = {'preimp': fpre,
             'postimp': fpost}
    temp = []
    for d, vcf in steps.items():
        df = pd.DataFrame.from_dict(vcfdf.PandasVCF(vcf, idt=idt).aaf,
                                    orient='index',
                                    columns=[d + '_' + set])
        temp.append(df)

    return temp


def get_df_evol(err_dic, fpre, fpost, path_all, gen, chk, idt='id'):
    bin_aaf = prm.BIN_AAF

    df_aaf = alltls.get_aaf(path_all, id=idt)

    for k_set in err_dic.keys():
        list_aaf = compute_aaf_evol(k_set, fpre, fpost, gtgl=gen, chk_sz=chk, idt=idt)
        df_aaf = df_aaf.join(list_aaf)
    df_aaf.sort_values(by='af_info', axis=0, inplace=True)

    if bin_aaf:
        convert = np.vectorize(lambda x: alltls.convert_aaf(x))
        inter = prm.INTER
        np_bin_aaf = np.digitize(convert(df_aaf.values),
                                 bins=inter)
        np_bin_aaf = np.subtract(np_bin_aaf/len(inter), 0.00)
        df_aaf = pd.DataFrame(data=np_bin_aaf, index=df_aaf.index, columns=df_aaf.columns)

    return df_aaf


def plot_test(df_aaf, col_set, typ='scatter'):
    plt.rcParams["figure.figsize"] = [12, 6]
    plt.rcParams["figure.autolayout"] = True

    lin_aaf, ax_lin = plt.subplots()
    for a, k_set in enumerate(col_set):
        colors = tuple(np.random.rand(3))  # RGBA --> 3
        if typ == 'scatter':
            df_aaf.plot.scatter(x='af_info',
                                y='preimp_' + k_set,
                                ax=ax_lin,
                                label='preimp_' + k_set,
                                marker='v',
                                color=colors)
            df_aaf.plot.scatter(x='af_info',
                                y='postimp_' + k_set,
                                ax=ax_lin,
                                label='postimp_' + k_set,
                                marker='o',
                                color=colors)
        if typ == 'line':
            df_aaf.plot(x='af_info',
                        y='preimp_' + k_set,
                        ax=ax_lin,
                        linestyle='--',
                        color=colors)
            df_aaf.plot(x='af_info',
                        y='postimp_' + k_set,
                        ax=ax_lin,
                        linestyle='-',
                        color=colors)
        plt.xlabel('Theoretical AAF')
        plt.ylabel('Dataset AAF')
        plt.plot(range(2), linestyle='-', color='k')
        delta_post = df_aaf['postimp_' + k_set].sub(df_aaf['af_info'])
        err_post = alltls.rmse_df(delta_post.to_frame(), kind='mse')
        plt.text(0.65, 0.25-0.1*a, 'MSE postimp_' + k_set + ' = ' + str(err_post))
        plt.legend()

    plt.title('AAF evolution through processing of data sets', loc='center')
    plt.suptitle("")
    plt.savefig('aaf_evol.{}.knockedout.pooled.chunk{}.png'.format(typ, df_aaf.shape[0]),
                orientation='landscape',
                dpi='figure')


def plot_1ko_vs_1batch(df_aaf, typ='scatter'):
    plt.rcParams["figure.figsize"] = [12, 6]
    plt.rcParams["figure.autolayout"] = True

    lin_aaf, ax_lin = plt.subplots()
    colors = tuple(np.random.rand(3))  # RGBA --> 3
    annot = zip(df_aaf.index,
                zip(df_aaf.loc[:, 'postimp_' + 'one_batch'],
                    df_aaf.loc[:, 'postimp_' + 'one_knocked_out']
                    )
                )

    if typ == 'scatter':
        # df_aaf.plot.scatter(x='preimp_' + 'one_batch',
        #                     y='preimp_' + 'one_knocked_out',
        #                     ax=ax_lin,
        #                     label='preimp_' + 'one_knocked_out',
        #                     marker='v',
        #                     color=colors)
        df_aaf.plot.scatter(x='postimp_' + 'one_batch',
                            y='postimp_' + 'one_knocked_out',
                            ax=ax_lin,
                            label='postimp_' + 'one_knocked_out',
                            marker='o',
                            color=colors)
        my_annotate2 = partial(ax_lin.annotate,
                               fontsize=6)
        all(starmap(my_annotate2, annot))  # all() for unpacking iterator over a function, NOT list
        aaf_leave = df_aaf.loc[:, 'postimp_one_knocked_out']
        aaf_batch = df_aaf.loc[:, 'postimp_one_batch']
        aaf_diff = aaf_batch.sub(aaf_leave).abs()
        df_aaf['aaf_diff'] = aaf_diff
        print('Outlier markers:\n', df_aaf.query('aaf_diff > 0.2').loc[:, ['af_info',
                                                                           'postimp_one_knocked_out',
                                                                           'postimp_one_batch',
                                                                           'aaf_diff']])

    if typ == 'line':
        # df_aaf.plot(x='preimp_' + 'one_batch',
        #             y='preimp_' + 'one_knocked_out',
        #             ax=ax_lin,
        #             linestyle='--',
        #             color=colors)
        df_aaf.plot(x='postimp_' + 'one_batch',
                    y='postimp_' + 'one_knocked_out',
                    ax=ax_lin,
                    linestyle='-',
                    color=colors)

    plt.xlabel('One-Batch AAF')
    plt.ylabel('One-Knocked-Out AAF')
    plt.plot(range(2), linestyle='-', color='k')
    plt.legend()
    plt.title('AAF evolution through processing of data sets', loc='center')
    plt.suptitle("")
    plt.savefig('aaf_evol.{}.1ko_vs_1batch.pooled.chunk{}.png'.format(typ, df_aaf.shape[0]),
                orientation='landscape',
                dpi='figure')



if __name__ == '__main__':
    if prm.GTGL == 'GT':
        cd = prm.PATH_GT_FILES
    if prm.GTGL == 'GL':
        if prm.unknown_gl != 'adaptative':
            cd = os.path.join(prm.WD, 'gl', 'gl_' + '_'.join(np.around(prm.unknown_gl, decimals=2).astype(str)))
        else:
            cd = os.path.join(prm.WD, 'gl', 'gl_adaptative')

    folders = {'one_knocked_out': 'ko_markers_merged',
               'one_batch': 'all_snps_all_samples'}

    print('Load parameters'.ljust(80, '.'))
    sorting = True  # sort data sets by AAF and population values
    os.chdir(cd)

    params = [('gt', 'gt', prm.CHK_SZ, 'chrom:pos'),
              ('gl', 'one_knocked_out', prm.CHK_SZ, 'chrom:pos'),
              ('gl', 'one_batch', prm.CHK_SZ, 'chrom:pos')
              ]
    devol = []

    # Load AAFs
    pdvcfall = vcfdf.PandasVCF(os.path.join(prm.PATH_GT_FILES,
                                             'ALL.chr20.pooled.snps.gt.chunk{}.vcf.gz'.format(prm.CHK_SZ),
                                             indextype='id'))
    allaafs = pdvcfall.concatcols([pdvcfall.af_info, pdvcfall.aaf])
    afinfo = pdvcfall.af_info.to_frame()
    pdvcfimp = vcfdf.PandasVCF(os.path.join(prm.PATH_GT_FILES,
                                             '/IMP.chr20.pooled.snps.gt.chunk{}.vcf.gz'.format(prm.CHK_SZ),
                                             idt='id'))
    impaafs = pdvcfimp.concatcols([pdvcfimp.af_info, pdvcfimp.aaf])
    compaafs = allaafs.join(impaafs, how='inner', rsuffix='_imp')

    for p_set in params:
        print('\nSet params:', p_set)
        gtgl, batch, sz, idt = p_set
        ALL = prm.PATH_GT_FILES + '/ALL.chr20.snps.gt.chunk{}.vcf.gz'.format(sz)
        B1 = prm.PATH_GT_FILES + '/IMP.chr20.pooled.snps.gt.chunk{}.vcf.gz'.format(sz)
        if gtgl == 'gl':
            POOL = os.path.join(cd,
                                folders[batch],
                                'IMP.chr20.pooled.imputed.gt.chunk{}.vcf.gz'.format(sz))
        if gtgl == 'gt':
            POOL = prm.PATH_GT_FILES + '/IMP.chr20.pooled.imputed.gt.chunk{}.vcf.gz'.format(sz)

        print('Raw samples and index labels'.ljust(80, '.'))
        # Create line-index and column-index (multi-indexing)
        aaf_idx, pop_idx = alltls.make_index(B1, src=ALL, id=idt)

        print('Build datasets'.ljust(80, '.'))
        rawvcf = vcfdf.PandasVCF(B1, indextype=idt)
        poolvcf = vcfdf.PandasVCF(POOL, indextype=idt)
        raw0, raw1 = rawvcf.vcf2dframe()
        pool0, pool1 = poolvcf.vcf2dframe()
        # raw0 = raw0.join(pool0.index.to_frame(), how='inner').iloc[:, :-1]
        # raw1 = raw1.join(pool1.index.to_frame(), how='inner').iloc[:, :-1]

        raw0 = raw0.join(afinfo, how='inner').drop('af_info', axis=1)
        raw1 = raw1.join(afinfo, how='inner').drop('af_info', axis=1)
        pool0 = pool0.join(afinfo, how='inner').drop('af_info', axis=1)
        pool1 = pool1.join(afinfo, how='inner').drop('af_info', axis=1)

        samples = VCF(B1).samples
        variants = raw0.index.tolist()

        # NumPy juxtapos for error computation
        raw = np.stack([raw0.values, raw1.values], axis=-1)
        pool = np.stack([pool0.values, pool1.values], axis=-1)

        set_errors = alltls.compute_imp_err(['pooled'],
                                            [pool],
                                            raw,
                                            aaf_idx,
                                            pop_idx,
                                            sorting)
        df_chk = get_df_evol(set_errors,
                             B1,
                             POOL,
                             ALL,
                             gtgl,
                             sz,
                             idt=idt).loc[:, ['preimp_pooled', 'postimp_pooled']]
        df_chk.rename(columns={'preimp_pooled': 'preimp_{}'.format(batch),
                               'postimp_pooled': 'postimp_{}'.format(batch)},
                      inplace=True)

        devol.append(df_chk)

    df_plot = afinfo.join(devol, how='inner')
    df_plot.sort_values(by='af_info', inplace=True)

    plot_test(df_plot, col_set=list(zip(*params))[1], typ='line')
    plot_test(df_plot, col_set=list(zip(*params))[1], typ='scatter')
    plot_1ko_vs_1batch(df_plot, typ='scatter')

    candidates = alltls.extract_variant_onfreq(os.path.join(prm.PATH_GT_FILES,
                                                            prm.RAW['ref'].replace('.gl', '.gt')),
                                               [0.495, 0.505])
    snp_to_rm: list = candidates['id'].str.strip(' ').str.slice(3).values
    df_ko_only = df_plot.loc[['20:' + snp for snp in snp_to_rm]]
    plot_test(df_ko_only, col_set=list(zip(*params))[1], typ='line')
    plot_test(df_ko_only, col_set=list(zip(*params))[1], typ='scatter')
    plot_1ko_vs_1batch(df_ko_only, typ='scatter')
"""
Outlier markers:
              
"""