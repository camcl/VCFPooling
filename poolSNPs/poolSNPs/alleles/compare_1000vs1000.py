import itertools
from scripts.poolSNPs.alleles import alleles_tools as alltls
from scripts.poolSNPs.results import *
from scripts.poolSNPs import parameters as prm
from persotools.files import *
from persotools.debugging import *

dbg = MyPrintClass(True)


def get_df_evol(err_dic, path_all, chk):
    bin_maf = prm.BIN_MAF

    df_maf = alltls.get_maf(path_all)

    for k_set in err_dic.keys():
        list_maf = alltls.compute_maf_evol(k_set, chk_sz=chk)
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


def plot_test(df_maf, col_set):
    plt.rcParams["figure.figsize"] = [12, 6]
    plt.rcParams["figure.autolayout"] = True

    lin_maf, ax_lin = plt.subplots()
    for a, k_set in enumerate(col_set):
        colors = tuple(np.random.rand(4))
        df_maf.plot.line(x='maf',
                         y='preimp_' + k_set,
                         ax = ax_lin,
                         linestyle='--',
                         color=colors)
        df_maf.plot.line(x='maf',
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
    plt.savefig('maf_evol.both.pooled.png',
                orientation='landscape',
                dpi='figure')


if __name__ == '__main__':
    ### DATA
    print('Load parameters'.ljust(80, '.'))
    sorting = True # sort data sets by MAF and population values
    os.chdir(prm.WD)
    params = ((True, 10000, 1000),
              (False, 1000, None))
    devol = []

    # Load MAFs
    mafs = alltls.get_maf('IMP.chr20.beagle1.chunk10000.vcf.gz')
    mafs = mafs.loc[:, 'maf'].to_frame()

    for p_set in params:
        print('\nSet params:', p_set)
        subset, sz, subchunk = p_set
        ALL = 'ALL.chr20.snps.gt.chunk{}.vcf.gz'.format(sz)
        B1 = 'IMP.chr20.beagle1.chunk{}.vcf.gz'.format(sz)
        POOL = 'IMP.chr20.pooled.beagle2.gt.chunk{}.corr.vcf.gz'.format(sz)

        print('Raw samples and index labels'.ljust(80, '.'))
        r_size, ms_size = len(VCF(B1).samples), len(VCF(POOL).samples)

        # Create line-index and column-index (multi-indexing)
        maf_idx, pop_idx = make_index(B1, src=ALL)

        ### PROCESSING
        print('Build datasets'.ljust(80, '.'))
        raw0, raw1 = alltls.vcf2dframe(VCF(B1), r_size)
        pool0, pool1 = alltls.vcf2dframe(VCF(POOL), ms_size)

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
                             ALL,
                             sz).loc[:, ['preimp_pooled', 'postimp_pooled']]
        df_chk.rename(columns={'preimp_pooled': 'preimp_' + str(sz),
                               'postimp_pooled': 'postimp_' + str(sz)},
                      inplace=True)

        devol.append(df_chk)

    df_plot = mafs.join(devol, how='inner')
    df_plot.sort_values(by='maf', inplace=True)
    print('dfplot\n', df_plot)
    plot_test(df_plot, col_set=['10000', '1000'])
