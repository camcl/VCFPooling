import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import seaborn as sns
from scipy.stats import *
import numpy as np
import pandas as pd
from . import alleles_tools as alltls

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
    im = ax.imshow(dtfr.values)

    ax.set_xticks(np.arange(dtfr.shape[1]))
    ax.set_yticks(np.arange(dtfr.shape[0]))
    ax.set_xticklabels(col)
    ax.set_yticklabels(idx)
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
        ay.set_yticks(np.arange(len(idY)))
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


def plot_mafs(err_dic, boxplot=True, lineplot=True):
    if boxplot:
        box_maf, ax_box = plt.subplots(nrows=1, ncols=3, sharey=True)
    if lineplot:
        lin_maf, ax_lin = plt.subplots()
    df_maf = alltls.get_maf('ALL.chr20.snps.gt.chunk.vcf.gz').set_index('id')
    a = 0
    colors = ['b', 'g', 'r']
    dic_maf = {}
    for k_set in err_dic.keys():
        dic_maf = alltls.compute_maf_evol(k_set)
        df_maf = pd.concat([df_maf,
                            pd.DataFrame.from_dict(dic_maf)],
                           axis=1,
                           sort=True)
        df_maf.sort_values(by='maf', axis=0, inplace=True)

        #print(df_maf.tail(10))
    #
    #     if boxplot:
    #         df_maf['delta_pre_' + k_set] = (df_maf['preimp_' + k_set] - df_maf['maf']).map(abs)
    #         df_maf['delta_post_' + k_set] = (df_maf['postimp_' + k_set] - df_maf['maf']).map(abs)
    #         df_maf[['delta_pre_' + k_set, 'delta_post_' + k_set]].plot.box(by='maf_inter',
    #                                                                    ax=ax_box[a],
    #                                                                    showfliers=False,
    #                                                                    showmeans=True)
    #
    #         plt.ylim(0.0, 0.001)
    #
    #     if lineplot:
    #         df_maf.plot.line(x='maf',
    #                          y='preimp_' + k_set,
    #                          ax = ax_lin,
    #                          linestyle='--',
    #                          color=colors[a])
    #         df_maf.plot.line(x='maf',
    #                             y='postimp_' + k_set,
    #                             ax=ax_lin,
    #                             linestyle='-',
    #                             color=colors[a])
    #         plt.plot(range(2), linestyle='-', color='k')
    #         plt.legend()
    #
    #     a += 1
    # plt.title('MAF evolution through processing of data sets', loc='center')
    # plt.suptitle("")
    # plt.show()