import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
import pandas as pd
from scripts.VCFPooling.python.archived.alleles import alleles_tools as alltls
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


def plot_aaf_twist(setname: str, setdf: pd.DataFrame, vcfpath: str):
    # TODO: Refactor: not any longer needed if overlay in plot_aaf_ scripts
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


def plot_aaf_vs_miss():
    #todo: transform into a separate script?
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