import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from cyvcf2 import VCF

from scripts.VCFPooling.python.archived.alleles import alleles_tools as alltls
from scripts.VCFPooling.poolSNPs import dataframe as vcfdf
from scripts.VCFPooling.poolSNPs import parameters as prm
from persotools.files import *

"""
Plot correlation between heterozygosity/homozygosity for every marker
 after imputation (computed with Beagle)
and heterozygosity/homozygosity in the original dataset (theoretical).
Compare effects of imputation:
- GT vs GL
- varying GL values for filling in missing genotypes
"""


def get_df_evol(fpre, fpost, dset, idt='chrom:pos'):
    #df_zyg = alltls.get_aaf(path_all, id=idt)

    for z in ['het', 'hom_alt', 'hom_ref']:
        list_zyg = alltls.compute_zygosity_evol(z, fpre, fpost, idt=idt)
        df_zyg = pd.DataFrame.join([dfz.add_suffix('_' + dset) for dfz in list_zyg])

    return df_zyg


def plot_zygosity_correlation(df_zyg, col_set, typ='scatter'):
    plt.rcParams["figure.figsize"] = [12, 6]
    plt.rcParams["figure.autolayout"] = True

    lin_zyg, ax_lin = plt.subplots()
    for a, z_set in enumerate(col_set):
        colors = tuple(numpy.random.rand(3)) # RGBA --> 3
        if typ == 'scatter':
            df_zyg.plot.scatter(x='af_info',
                                y='preimp_' + z_set,
                                ax=ax_lin,
                                label='preimp_' + z_set,
                                marker='v',
                                color=colors)
            df_zyg.plot.scatter(x='af_info',
                                y='postimp_' + z_set,
                                ax=ax_lin,
                                label='postimp_' + z_set,
                                marker='o',
                                color=colors)
        if typ == 'line':
            df_zyg.plot(x='af_info',
                        y='preimp_' + z_set,
                        ax=ax_lin,
                        linestyle='--',
                        color=colors)
            df_zyg.plot(x='af_info',
                        y='postimp_' + z_set,
                        ax=ax_lin,
                        linestyle='-',
                        color=colors)
        plt.xlabel('Theoretical zyg')
        plt.ylabel('Dataset zyg')
        plt.plot(range(2), linestyle='-', color='k')
        delta_post = df_zyg['postimp_' + z_set].sub(df_zyg['preimp_' + z_set])
        err_post = alltls.rmse_df(delta_post.to_frame(), kind='mse')
        plt.text(0.65, 0.25-0.1*a, 'MSE postimp_' + z_set + ' = ' + str(err_post))
        plt.legend()

    plt.title('zyg evolution through processing of data sets', loc='center')
    plt.suptitle("")
    plt.savefig('zyg_evol.{}.gtgl.pooled.chunk{}.png'.format(typ, df_zyg.shape[0]),
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

    df_plot = afinfo.join(devol, how='inner')
    print(df_plot)
    df_plot.sort_values(by='af_info', inplace=True)
    print('dfplot\n', df_plot)
    plot_zygosity_correlation(df_plot, col_set=list(zip(*params))[0], typ='scatter')

