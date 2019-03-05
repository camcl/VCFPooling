import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import seaborn as sns
from scipy.stats import *
import math
import numpy as np
import pandas as pd
from cyvcf2 import VCF
from . import alleles_tools as alltls
from scripts.poolSNPs import parameters as prm
from persotools.debugging import *
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

    plt.rcParams["figure.autolayout"] =  True
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

    plt.savefig(figname +'.heatmap{}.png'.format('.sorted' if sorting else ''),
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
        lev_err.drop(['id', 'maf'], axis=1, inplace=True)
        lev_err.boxplot(by='maf_inter',
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
    plt.savefig('root.mean.square.error.box.density.png', dpi='figure')


def plot_mafs(err_dic, path_all, boxplot=True, lineplot=True):
    if boxplot:
        box_maf, ax_box = plt.subplots(nrows=1, ncols=3, sharey=True)
    if lineplot:
        lin_maf, ax_lin = plt.subplots()
    df_maf = alltls.get_maf(path_all).set_index('id')
    a = 0
    colors = ['b', 'g', 'r']
    for k_set in err_dic.keys():
        dic_maf = alltls.compute_maf_evol(k_set)
        df_maf = pd.concat([df_maf,
                            pd.DataFrame.from_dict(dic_maf)],
                           axis=1,
                           sort=True)
        df_maf.sort_values(by='maf', axis=0, inplace=True)

        if boxplot:
            df_maf['delta_pre_' + k_set] = (df_maf['preimp_' + k_set] - df_maf['maf']).map(abs)
            df_maf['delta_post_' + k_set] = (df_maf['postimp_' + k_set] - df_maf['maf']).map(abs)
            df_maf[['delta_pre_' + k_set, 'delta_post_' + k_set]].plot.box(by='maf_inter',
                                                                       ax=ax_box[a],
                                                                       showfliers=False,
                                                                       showmeans=True)

            plt.ylim(0.0, 0.001)

        if lineplot:
            df_maf.plot.line(x='maf',
                             y='preimp_' + k_set,
                             ax = ax_lin,
                             linestyle='--',
                             color=colors[a])
            df_maf.plot.line(x='maf',
                                y='postimp_' + k_set,
                                ax=ax_lin,
                                linestyle='-',
                                color=colors[a])
            plt.xlabel('Theoretical MAF')
            plt.ylabel('Dataset MAF')
            plt.plot(range(2), linestyle='-', color='k')
            plt.legend()

        a += 1
    plt.title('MAF evolution through processing of data sets', loc='center')
    plt.suptitle("")
    plt.savefig('maf_evol.png',
                orientation='landscape',
                dpi='figure')


def plot_box_low_maf(df_err, name,  ax):
    # cast maf_inter as a column
    cols = df_err.columns.droplevel(level='Population')
    df = pd.DataFrame(data=df_err.values, index=df_err.index, columns=cols)
    df.drop('maf', axis=1, inplace=True)
    rmse = alltls.rmse_df(df, ax=1).to_frame()
    rmse.reset_index(level='maf_inter', drop=False, inplace=True)

    ax = rmse.boxplot(by='maf_inter', ax=ax, showmeans=True, showfliers=False)
    plt.suptitle('RMSE for the {} data set'.format(name))
    return ax


def plot_err_vs_het(dset, err_set, file_in, err_kind='rmse', low_maf=False, save=True, ax=None):

    het = alltls.per_site_heteroz(file_in)
    ogn = 'imputed' if file_in.find('beagle') is not -1 else 'true'
    nb_samples = len(VCF(file_in).samples)

    if ax is None:
        fig, ax = plt.subplots()

    # Create rmse DataFrame with MultiIndex
    s = alltls.rmse_df(err_set, kind=err_kind, ax=1) # ax=1: apply along columns
    df = pd.DataFrame(data=s.values, index=err_set.index, columns=[err_kind])
    df.reset_index(level=['maf', 'maf_inter'], drop=False, inplace=True)

    # heterozygosity in the final input file
    htzSer = pd.Series(het[:,-1].astype(float), index=het[:, 0])
    htzPct = htzSer.apply(lambda h: h * 100 / nb_samples).rename('het')
    df['het'] = htzPct

    if not low_maf:
        ax = df.plot.scatter('het', err_kind, c='maf', cmap=plt.cm.viridis, ax=ax)
    else:
        df_low = df.query('maf_inter <= 2')
        ax = df_low.plot.scatter('het', err_kind, c='maf', cmap=plt.cm.PuBu, ax=ax)
    ax.set_title('{} vs. {} heterozygosity in {} dataset'.format(err_kind.upper(), ogn, dset))

    if save:
        ax = None
        # plt.axhline(y=site_err,
        #             linestyle='dashed',
        #             color='k',
        #             label='mean error')
        plt.savefig('{}.heterozygosity.{}.png'.format(err_kind, dset),
                    dpi='figure')

    return ax


def plot_err_vs_miss(k, err_set, err_kind='rmse'):

    nb_samples = err_set.shape[1]
    site_err = alltls.rmse_df(err_set, kind=err_kind, ax=1).rename(err_kind).to_frame()
    site_err.reindex(index=err_set.index, copy=False)
    site_err.reset_index(level=['maf'], drop=False, inplace=True)
    site_err['maf'].astype(float, copy=False)
    # missing data from the preimputed dataset
    misNp = alltls.count_missing_alleles('IMP.chr20.{}.cfgt.vcf.gz'.format(k))
    misSer = pd.Series(misNp[:, -1], index=misNp[:, 0]).astype(float, copy=True)
    misPct = misSer.apply(lambda h: h * 100 / nb_samples).rename('miss').to_frame()
    site_err = site_err.join(misPct, on='id')
    site_err.plot.scatter('miss', err_kind, c='maf', cmap=plt.cm.viridis)
    plt.title('Imputation error vs. missing data rate in {} data set'.format(k))
    plt.savefig('rmse.missing.data.{}.png'.format(k.upper()),
                dpi='figure')


def multiplot_err_het(err):

    plt.rcParams["figure.figsize"] = [12, 6]
    low_maf, axm = plt.subplots(3, 2)
    j = 0
    for k in ['missing', 'pooled']:
        i = 0
        ek = err[k].reset_index(level=['maf'], drop=False, inplace=False)
        low_err = ek.query('maf_inter < 3', inplace=False)
        print('Marker with maximum error: ', low_err.mean(axis=1).idxmax(axis=0))
        print(low_err.query('maf_inter == 2', inplace=False).mean(axis=1).describe())
        print(low_err.query('maf_inter == 1', inplace=False).mean(axis=1).describe())
        low_err.set_index('maf', drop=True, append=True, inplace=True)
        file1 = 'IMP.chr20.{}.beagle2.corr.vcf.gz'.format(k)
        axm[i,j] = plot_err_vs_het(k,
                                   low_err,
                                   file1,
                                   err_kind='mse',
                                   low_maf=True,
                                   save=False,
                                   ax=axm[i,j])
        i += 1
        file2 = prm.RAW['imp']
        axm[i, j] = plot_err_vs_het(k,
                                    low_err,
                                    file2,
                                    err_kind='rmse',
                                    low_maf=True,
                                    save=False,
                                    ax=axm[i,j])
        i += 1
        axm[i, j] = plot_err_vs_het(k,
                                    low_err,
                                    file2,
                                    err_kind='mse',
                                    low_maf=True,
                                    save=False,
                                    ax=axm[i,j])
        j += 1
    #plt.show()
    plt.savefig('root.mean.square.error.maf.low.png', dpi='figure')
