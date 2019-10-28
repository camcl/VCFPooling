import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.stats import *
import numpy as np
import pandas as pd
from cyvcf2 import VCF
from scripts.VCFPooling.poolSNPs.alleles import alleles_tools as alltls
from scripts.VCFPooling.poolSNPs import parameters as prm
from persotools.files import *


"""
Tools for plotting alleles and genotype data.
"""


def plot_heat_map(dtfr, figname, figsize, sorting, title='{}', rightYaxis=False):
    """
    Plots the values of a dataframe as a heatmap and save the figure as a PNG file.
    :param dtfr: dataframe with data to plot
    :param figname: string for the PNG name to save
    :param figsize: pair of integer values (list-like) giving width and height of the figure
    :param title: string for titling the figure
    :param rightYaxis: boolean for plotting a second y-axis or not
    :return:
    """
    #TODO: colorbar's parameters to tune
    plt.rcParams["figure.figsize"] = figsize

    try:
        idx = np.array(list(zip(*dtfr.index)))[0]
        idY = np.array(list(zip(*dtfr.index)))[-1]
        col = np.array(list(zip(*dtfr.columns)))[0]

    except:
        idx = dtfr.index
        col = dtfr.columns

    plt.rcParams["figure.autolayout"] = True
    fig, ax = plt.subplots()
    im = ax.imshow(dtfr.transpose().values)

    ax.set_yticks(np.arange(dtfr.shape[1]))
    ax.set_xticks(np.arange(dtfr.shape[0]))
    ax.set_yticklabels(col)
    ax.set_xticklabels(idx)
    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(),
             fontsize=0.25,
             rotation=45,
             ha="center",
             rotation_mode="anchor")
    plt.setp(ax.get_yticklabels(),
             fontsize=0.25,
             va="center")

    if rightYaxis:
        ay = fig.add_subplot(111, sharex=ax, frameon=False)
        ay.set_yticks(np.arange(len(idx)))
        #ay.set_yticklabels(np.arange(1, max(idY) + 1).astype(str))
        ay.set_yticklabels(['1', '2', '3'])
        ay.yaxis.grid(which="major", color='w', linestyle='-', linewidth=0.2)
        ay.yaxis.tick_right()
        ay.yaxis.set_label_position("right")
        plt.setp(ay.get_yticklabels(),
                 fontsize=0.25,
                 va="center")

    ax.set_title(title.format(figname))
    axins = inset_axes(ax,
                        width=2,
                        height=0.1,
                        loc='lower center')
    plt.colorbar(im,
                 cax=axins,
                 orientation='horizontal',
                 fraction=0.01,
                 anchor=(0.5,-1.0))
                        # panchor=(0.5, 0.0))
    # cbar.ax.set_xlabel(legend,
    #                    va="bottom")
                       # fontsize=0.1)
    #fig.tight_layout(h_pad=0.5)

    plt.savefig(figname + '.heatmap{}.chunk{}.png'.format('.sorted' if sorting else '', prm.CHK_SZ),
                orientation='landscape',
                dpi='figure')


def boxplot_densities(set_errors):
    pd.set_option('precision', 25)
    plt.rcParams["figure.figsize"] = [20, 20]
    plt.rcParams["figure.autolayout"] = True
    errbox, axs = plt.subplots(3, 2)

    import sys
    sys.stdout = open('.log', 'w')

    pos = 0
    for k in ['pooled', 'missing']:
        df = set_errors[k]['grid']
        pop_err = df.groupby(axis=1, level='Population').mean()
        lev_err = alltls.rmse_df(df)
        lev_err.drop(['id', 'af_info'], axis=1, inplace=True)
        lev_err.boxplot(by='aaf_bin',
                        ax=axs[pos, 0],
                        showmeans=True,
                        showfliers=False)
        plt.suptitle("")
        axs[pos, 0].set_title('rmse for the {} data set'.format(k))

        print('\r\nNumber of missing genotypes in the "{}" dataset:'.format(k))
        print('Mean Error: {}'.format(set_errors[k]['mean']))

        axs[pos, 1].set_xlim(left=-1.0, right=1.0)
        pop_err.plot.density(ax=axs[pos, 1],
                             legend=False)
        pop_skw = skew(pop_err.values, axis=0, bias=True)
        axs[pos, 1].set_title('Error density distributions for the {} data set'.format(k))
        pos += 1
    handles, labels = axs[1, 1].get_legend_handles_labels()
    detailed_labels = [pop + ', skewness = ' + str(sk) for pop, sk in zip(labels, pop_skw)]
    axs[1, 1].legend(handles,
                     detailed_labels,
                     loc='center left',
                     bbox_to_anchor=(1.05, 0.5),
                     borderaxespad=0,
                     frameon=True)
    plt.savefig('root.mean.square.error.box.density.chunk{}.png'.format(prm.CHK_SZ), dpi='figure')


def plot_aaf_evol(err_dic, path_all, typ='line'):
    bin_aaf = prm.BIN_AAF

    plt.rcParams["figure.figsize"] = [12, 6]
    plt.rcParams["figure.autolayout"] = True

    pdvcf = alltls.PandasVCF(path_all, indextype='chrom:pos')
    df_aaf = pdvcf.concatcols([pdvcf.af_info, pdvcf.aaf])

    if prm.SUBSET:
        df_aaf = df_aaf.iloc[:prm.SUBCHUNK]

    for k_set in err_dic.keys():
        list_aaf = alltls.compute_aaf_evol(k_set)
        df_aaf = df_aaf.join(list_aaf)

    df_aaf.sort_values(by='af_info', axis=0, inplace=True)

    if bin_aaf:
        df_aaf = df_aaf.join([pdvcf.aaf_binned(b=prm.INTER)])

    lin_aaf, ax_lin = plt.subplots()
    a = 0
    colors = ['b', 'g', 'r']
    for k_set in err_dic.keys():
        if typ == 'line':
            df_aaf.plot.line(x='af_info',
                             y='preimp_' + k_set,
                             ax = ax_lin,
                             linestyle='--',
                             color=colors[a])
            df_aaf.plot.line(x='af_info',
                             y='postimp_' + k_set,
                             ax=ax_lin,
                             linestyle='-',
                             color=colors[a])
        if typ == 'scatter':
            df_aaf.plot.line(x='af_info',
                             y='preimp_' + k_set,
                             ax=ax_lin,
                             marker='v',
                             color=colors[a])
            df_aaf.plot.line(x='af_info',
                             y='postimp_' + k_set,
                             ax=ax_lin,
                             marker='o',
                             color=colors[a])
        plt.xlabel('Theoretical AAF')
        plt.ylabel('Dataset AAF')
        plt.plot(range(2), linestyle='-', color='k')
        delta_post = df_aaf['postimp_' + k_set].sub(df_aaf['af_info'])
        err_post = alltls.rmse_df(delta_post.to_frame(), kind='mse')
        plt.text(0.65, 0.25-0.1*a, 'MSE postimp_' + k_set + ' = ' + str(err_post))
        plt.legend()
        a += 1

    plt.title('AAF evolution through processing of data sets', loc='center')
    plt.suptitle("")
    plt.savefig('af_info_evol.{}.chunk{}.png'.format(typ, prm.CHK_SZ),
                orientation='landscape',
                dpi='figure')


def plot_aaf_twist(setname: str, setdf: pd.DataFrame, vcfpath: str):
    """
    Plot aaf before/after pooling.
    Plot error after pooling vs. aaf
    :param setname: pooled/missing
    :param setdf: dataframe with discordance after imputation
    :param vcfpath: path to vcf for aaf
    :param path_all:
    :return:
    """
    plt.rcParams["figure.figsize"] = [10, 6]
    plt.rcParams["figure.autolayout"] = True

    chkfile = os.path.join(prm.PATH_GT_FILES, prm.CHKFILE)
    vcfchk = alltls.PandasVCF(chkfile, indextype='id')
    chkinfo = vcfchk.af_info()

    # df_aaf = alltls.get_aaf(chkfile, idt='id')
    # df_aaf.drop(['aaf_bin'], axis=1, inplace=True)
    # df_aaf.drop(['aaf'], axis=1, inplace=True)

    # aafs = pd.DataFrame.from_dict(alltls.compute_aaf(vcfpath, idt='id'),
    #                               orient='index',
    #                               columns=['aaf'])
    aafs = alltls.PandasVCF(vcfpath, indextype='id').aaf()
    print(aafs)
    df_aaf = pd.concat([chkinfo, aafs], axis=1)

    imperrors = pd.DataFrame(setdf.mean(axis=1), columns=['error_' + setname])
    print('\r\nMean squared imputation error from {} data set = {}'.format(setname, imperrors.mean()))
    df_aaf = df_aaf.join(imperrors)
    # print(df_aaf)

    df_aaf.sort_values(by='af_info', axis=0, inplace=True)
    bins = np.linspace(0.0, 1.0, 10)
    x_pos = np.add(bins, 0.05).round(2)
    groups = pd.Series(np.digitize(df_aaf['af_info'].astype(float), bins=bins))
    df_aaf = df_aaf.assign(aaf_10bins=lambda x: np.digitize(df_aaf['af_info'].astype(float), bins=bins))
    df_aaf = df_aaf.assign(af_info10=lambda x: np.multiply(df_aaf['af_info'].astype(float), 10))

    a = 0
    colors = [['b', 'tab:blue'], ['g', 'tab:green'], ['r', 'tab:red']]
    print('Plotting for the {} data set'.format(setname).ljust(80, '.'))
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    df_aaf.plot.scatter(x='af_info',
                        y='aaf',
                        ax=ax1,
                        marker='o',
                        s=0.7,
                        color=colors[a][0],
                        label='aaf_{}_imputed'.format(setname))
    df_aaf.boxplot(column='error_' + setname,
                   by='aaf_10bins',
                   ax=ax2,
                   widths=0.05,
                   positions=x_pos,
                   showfliers=False,
                   grid=False
                   )
    ax1.plot([0.0, 1.0], [0.0, 1.0], linestyle='--', color='k', label='identity function')

    ax1.set_xlabel('Theoretical AAF')
    ax1.set_ylabel('Dataset AAF')
    ax1.set_xbound(lower=0.0, upper=1.0)
    ax2.set_ylabel('Mean imputation error (discordance rate)')
    ax2.set_title('')
    ax2.set_xbound(lower=0.0, upper=1.0)
    ax2.get_xaxis().set_ticklabels(x_pos)
    fig.suptitle('')
    plt.savefig('aaf_twist_imputation_error.{}.chunk{}.png'.format(setname, prm.CHK_SZ))
    plt.show()


def plot_err_vs_het(dset, err_set, file_in, err_kind='rmse', low_aaf=False, save=True, ax=None):

    het = alltls.per_site_heteroz(file_in)
    ogn = 'imputed' if file_in.find('beagle') is not -1 else 'true'
    nb_samples = len(VCF(file_in).samples)

    if ax is None:
        fig, ax = plt.subplots()

    # Create rmse DataFrame with MultiIndex
    s = alltls.rmse_df(err_set, kind=err_kind, ax=1) # ax=1: apply along columns
    df = pd.DataFrame(data=s.values, index=err_set.index, columns=[err_kind])
    df.reset_index(level=['af_info', 'aaf_bin'], drop=False, inplace=True)

    # heterozygosity in the final input file
    htzSer = pd.Series(het[:,-1].astype(float), index=het[:, 0])
    htzPct = htzSer.apply(lambda h: h * 100 / nb_samples).rename('het')
    df['het'] = htzPct

    if not low_aaf:
        ax = df.plot.scatter('het', err_kind, c='af_info', cmap=matplotlib.cm.viridis, ax=ax)
    else:
        df_low = df.query('aaf_bin <= 2')
        ax = df_low.plot.scatter('het', err_kind, c='af_info', cmap=matplotlib.cm.PuBu, ax=ax)
    ax.set_title('{} vs. {} heterozygosity in {} dataset'.format(err_kind.upper(), ogn, dset))

    if save:
        ax = None
        # plt.axhline(y=site_err,
        #             linestyle='dashed',
        #             color='k',
        #             label='mean error')
        plt.savefig('{}.heterozygosity.{}.chunk{}.png'.format(err_kind, dset, prm.CHK_SZ),
                    dpi='figure')

    return ax


def multiplot_err_het(err):

    plt.rcParams["figure.figsize"] = [12, 6]
    low_aaf, axm = plt.subplots(3, 2)
    j = 0
    for k in ['missing', 'pooled']:
        i = 0
        ek = err[k].reset_index(level=['af_info'], drop=False, inplace=False)
        low_err = ek.query('aaf_bin < 3', inplace=False)
        # print('Marker with maximum error: ', low_err.mean(axis=1).idxmax(axis=0))
        # print(low_err.query('aaf_bin == 2', inplace=False).mean(axis=1).describe())
        # print(low_err.query('aaf_bin == 1', inplace=False).mean(axis=1).describe())
        low_err.set_index('af_info', drop=True, append=True, inplace=True)
        file1 = (prm.POOLED['b2'] if k == 'pooled' else prm.MISSING['b2']) + '.vcf.gz'
        axm[i,j] = plot_err_vs_het(k,
                                   low_err,
                                   file1,
                                   err_kind='mse',
                                   low_aaf=True,
                                   save=False,
                                   ax=axm[i,j])
        i += 1
        file2 = prm.RAW['imp']
        axm[i, j] = plot_err_vs_het(k,
                                    low_err,
                                    file2,
                                    err_kind='rmse',
                                    low_aaf=True,
                                    save=False,
                                    ax=axm[i,j])
        i += 1
        axm[i, j] = plot_err_vs_het(k,
                                    low_err,
                                    file2,
                                    err_kind='mse',
                                    low_aaf=True,
                                    save=False,
                                    ax=axm[i,j])
        j += 1
    plt.savefig('root.mean.square.error.aaf.low.chunk{}.png'.format(prm.CHK_SZ), dpi='figure')


def plot_err_vs_miss(k, err_set, err_kind='rmse'):

    nb_samples = err_set.shape[1]
    site_err = alltls.rmse_df(err_set, kind=err_kind, ax=1).rename(err_kind).to_frame()
    site_err.reindex(index=err_set.index, copy=False)
    site_err.reset_index(level=['af_info'], drop=False, inplace=True)
    site_err['af_info'].astype(float, copy=False)

    # missing data from the preimputed dataset
    misNp = alltls.count_missing_alleles('IMP.chr20.{}.snps.chunk{}.vcf.gz'.format(k, prm.CHK_SZ))
    misSer = pd.Series(misNp[:, -1], index=misNp[:, 0]).astype(float, copy=True)
    misPct = misSer.apply(lambda h: h * 100 / nb_samples).rename('miss').to_frame()

    site_err = site_err.join(misPct, on='id')
    site_err.plot.scatter('miss', err_kind, c='af_info', cmap=matplotlib.cm.viridis)
    plt.title('Imputation error vs. missing data rate in {} data set'.format(k))
    plt.savefig('rmse.missing.data.{}.chunk{}.png'.format(k, prm.CHK_SZ),
                dpi='figure')


def plot_aaf_vs_miss():
    """
    Relevant for GT only
    :param k:
    :param err_set:
    :return:
    """
    if prm.GTGL == 'GT':
        os.chdir(prm.PATH_GT_FILES)
        chkgtfile = 'IMP.chr20.snps.gt.chunk{}.vcf.gz'.format(prm.CHK_SZ)
        # missing data from the preimputed dataset
        aafs = alltls.PandasVCF(chkgtfile, indextype='id').af_info()
        df_plot = aafs.to_frame()
        for (dic, name) in [(prm.POOLED, 'pooled'), (prm.MISSING, 'missing')]:
            misNp = alltls.count_missing_alleles(dic['imp'], id='id')
            misSer = pd.Series(misNp[:, -1], index=misNp[:, 0]).astype(float, copy=True)
            misPct = misSer.apply(lambda h: h * 100 / prm.NB_IMP).rename('miss_' + name).to_frame()
            df_plot = df_plot.join(misPct, how='inner')
        print(df_plot)

        fig, axis = plt.subplots()
        df_plot.plot.scatter(x='af_info',
                             y='miss_pooled',
                             s=2,
                             ax=axis,
                             label='pooled dataset',
                             marker='o',
                             color='tab:olive')
        df_plot.plot.scatter(x='af_info',
                             y='miss_missing',
                             s=2,
                             ax=axis,
                             label='random-deleted dataset',
                             marker='o',
                             color='tab:brown')
        axis.legend(loc='center', fontsize=6, bbox_to_anchor=(0.5, -0.3))
        plt.xlabel('AAF from VCF AF INFO-field')
        plt.ylabel('Data missing rate')
        fig.tight_layout()
        plt .savefig(prm.PLOTS_PATH + '/plot_aaf_vs_miss.chunk{}.jpg'.format(prm.CHK_SZ),
                     dpi=500)

    else:
        print('Function should be run with prm.GTGL = GT')


def plot_aaf_gl():
    """
    Shows how AAF = f impacts GL values for GT RR|RA|AA.
    :return:
    """
    import matplotlib.pyplot as plt

    f = np.arange(0.0, 1.0, 0.01)
    rr = np.vectorize(lambda x: (1 - x)**2)
    ra = np.vectorize(lambda x: 2 * (1 - x) * x)
    aa = np.vectorize(lambda x: x ** 2)

    plt.plot(f, rr(f), 'b-', label='RR')
    plt.plot(f, ra(f), 'g-', label='RA')
    plt.plot(f, aa(f), 'r-', label='AA')
    plt.legend()
    plt.show()


if __name__ == '__main__':
    plot_aaf_vs_miss()