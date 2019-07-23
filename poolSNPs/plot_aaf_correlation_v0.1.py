import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from cyvcf2 import VCF

from scripts.poolSNPs.alleles import alleles_tools as alltls
import scripts.poolSNPs.results as res
from scripts.poolSNPs import parameters as prm
from persotools.files import *

"""
Plot correlation between AAF after imputation (computed with Beagle)
and AAF in the original dataset (theoretical).
Compare effects of imputation:
- GT vs GL
- varying GL values for filling in missing genotypes
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
        df = pd.DataFrame.from_dict(alltls.compute_aaf(vcf, idt=idt),
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
        colors = tuple(np.random.rand(3)) # RGBA --> 3
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
    plt.savefig('aaf_evol.{}.gtgl.pooled.chunk{}.png'.format(typ, df_aaf.shape[0]),
                orientation='landscape',
                dpi='figure')


if __name__ == '__main__':
    ### DATA
    if prm.GTGL == 'GT':
        cd = prm.PATH_GT_FILES
    if prm.GTGL == 'GL':
        if prm.unknown_gl != 'adaptative':
            cd = os.path.join(prm.WD, 'gl', 'gl_' + '_'.join(np.around(prm.unknown_gl, decimals=2).astype(str)))
        else:
            # cd = os.path.join(prm.WD, 'gl')
            cd = os.path.join(prm.WD, 'gl', 'gl_adaptative', 'all_snps_all_samples')
    print('Load parameters'.ljust(80, '.'))
    sorting = True # sort data sets by AAF and population values
    os.chdir(cd)
    params = [('gt', prm.CHK_SZ, 'chrom:pos'), ('gl', prm.CHK_SZ, 'chrom:pos')]
    devol = []

    # Load AAFs
    allaafs = alltls.get_aaf(prm.PATH_GT_FILES + '/ALL.chr20.pooled.snps.gt.chunk{}.vcf.gz'.format(prm.CHK_SZ),
                             id='chrom:pos')
    aafs = allaafs.loc[:, 'af_info'].to_frame()
    print(aafs.shape)
    # aafs = allaafs.loc[:, 'aaf'].to_frame()
    impaafs = alltls.get_aaf(prm.PATH_GT_FILES + '/IMP.chr20.pooled.snps.gt.chunk{}.vcf.gz'.format(prm.CHK_SZ),
                             id='chrom:pos')
    compaafs = allaafs.loc[:, ['af_info', 'aaf']].join(impaafs.loc[:, ['af_info', 'aaf']], how='inner', rsuffix='_imp')

    for p_set in params:
        print('\nSet params:', p_set)
        gtgl, sz, idt = p_set
        ALL = prm.PATH_GT_FILES + '/ALL.chr20.snps.gt.chunk{}.vcf.gz'.format(sz)
        B1 = prm.PATH_GT_FILES + '/IMP.chr20.pooled.snps.gt.chunk{}.vcf.gz'.format(sz)
        if gtgl == 'gl':
            POOL = 'IMP.chr20.pooled.imputed.gt.chunk{}.vcf.gz'.format(sz)
        if gtgl == 'gt':
            POOL = os.path.join(prm.PATH_GT_FILES,
                                'all_snps_all_samples',
                                'IMP.chr20.pooled.imputed.gt.chunk{}.vcf.gz'.format(sz))

        print('Raw samples and index labels'.ljust(80, '.'))
        # Create line-index and column-index (multi-indexing)
        aaf_idx, pop_idx = alltls.make_index(B1, src=ALL, idt=idt)

        ### PROCESSING
        print('Build datasets'.ljust(80, '.'))
        raw0, raw1 = alltls.vcf2dframe(VCF(B1), prm.NB_IMP, id=idt)
        pool0, pool1 = alltls.vcf2dframe(VCF(POOL), prm.NB_IMP, id=idt)
        # raw0 = raw0.join(pool0.index.to_frame(), how='inner').iloc[:, :-1]
        # raw1 = raw1.join(pool1.index.to_frame(), how='inner').iloc[:, :-1]

        raw0 = raw0.join(aafs, how='inner').drop('af_info', axis=1)
        raw1 = raw1.join(aafs, how='inner').drop('af_info', axis=1)
        pool0 = pool0.join(aafs, how='inner').drop('af_info', axis=1)
        pool1 = pool1.join(aafs, how='inner').drop('af_info', axis=1)

        samples = VCF(B1).samples
        variants = raw0.index.tolist()

        # NumPy juxtapos for error computation
        raw = np.stack([raw0.values, raw1.values], axis=-1)
        pool = np.stack([pool0.values, pool1.values], axis=-1)
        print('shape raw set: ', raw.shape)
        print('shape pooled set: ', pool.shape)

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
        df_chk.rename(columns={'preimp_pooled': 'preimp_' + str(gtgl),
                               'postimp_pooled': 'postimp_' + str(gtgl)},
                      inplace=True)

        devol.append(df_chk)

    df_plot = aafs.join(devol, how='inner')
    print(df_plot)
    df_plot.sort_values(by='af_info', inplace=True)
    print('dfplot\n', df_plot)
    plot_test(df_plot, col_set=list(zip(*params))[0], typ='line')
    plot_test(df_plot, col_set=list(zip(*params))[0], typ='scatter')

