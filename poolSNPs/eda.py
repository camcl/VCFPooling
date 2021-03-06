import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
import sys, os
import argparse
import numpy as np
from scipy import stats
import pandas as pd
import seaborn as sns

# force PYTHONPATH to look into the project directory for modules
rootdir = os.path.dirname(os.path.dirname(os.getcwd()))
sys.path.insert(0, rootdir)

from VCFPooling.poolSNPs import dataframe as vcfdf

"""
Plot Early Data Analysis from pooling simulation:
The baseline is OmniExpress LD STU population + HD ref PAN (missing data, genotypes distribution)
vs.
OmniExpress HD in pooled STU population + HD ref PAN
- hypergeometric distortion
- missing data rate
- genotypes het/hom_ref/hom_alt proportions

To be included in BMC Bioinformatics article

Command line usage (assuming the current directory is VCFPooling/examples)
$ python3 -u ../poolSNPs/eda.py TEST.chr20.snps.gt.vcf.gz TEST.chr20.pooled.snps.gt.vcf.gz
$ python3 -u ../poolSNPs/eda.py IMP.chr20.snps.gt.vcf.gz IMP.chr20.pooled.snps.gt.vcf.gz
$ python3 -u ../poolSNPs/eda.py /home/camille/PoolImpHuman/data/20200812/IMP.chr20.missingHD.fullLD.snps.gt.vcf.gz /home/camille/PoolImpHuman/data/20200812/IMP.chr20.pooled.snps.gt.vcf.gz
"""

### COMMAND-LINE PARSING AND PARAMETERS
parser = argparse.ArgumentParser(description='Early Data Analysis with plots'
                                             'for the genetic structure  in a population'
                                             'of simulated pooled genotypes')
parser.add_argument('truegenos', metavar='trug', type=str, help='File with true genotypes GT', default=None)
parser.add_argument('pooledgenos', metavar='poog', type=str, help='File with pooled genotypes decoded into GT', default=None)
argsin = parser.parse_args()

paths = {'gt': {
    'true': argsin.truegenos,
    'pooled': argsin.pooledgenos}
}

dftrue = vcfdf.PandasMixedVCF(paths['gt']['true'], format='GT')
dfpool = vcfdf.PandasMixedVCF(paths['gt']['pooled'], format='GT')

df1 = dftrue.concatcols([dftrue.af_info, dftrue.missing_rate, dftrue.aaf])
df2 = dfpool.concatcols([dfpool.af_info, dfpool.missing_rate, dfpool.aaf])
df3 = dftrue.concatcols([dftrue.af_info, dftrue.het_rate, dftrue.hom_alt_rate, dftrue.hom_ref_rate])
df4 = dfpool.concatcols([dfpool.af_info, dfpool.het_rate, dfpool.hom_alt_rate, dfpool.hom_ref_rate])

print(df1.head(10))
print(df2.head(10))


# Plots' layout params

pts_sz = 10  # scatter dots size
pts_alpha = 1.0  # transparency
figsize = 8
labfont = 10
x_bins = [0.0, 0.02, 0.04, 0.06, 0.1, 0.2, 0.4, 0.6, 0.8, 0.9, 0.94, 0.96, 0.98, 1.0]
lab_bins = [0.01, 0.03, 0.05, 0.08, 0.15, 0.3, 0.5, 0.7, 0.85, 0.92, 0.95, 0.97, 0.99]
spaced2_labs = [0.01, None, 0.05, None, 0.15, None, 0.5, None, 0.85, None, 0.95, None, 0.99]
spaced3_labs = [None, 0.03, None, 0.08, None, 0.3, None, 0.7, None, 0.92, None, 0.97, None]
plt.rc('legend', fontsize=7)  # legend texts size


# Functions


def inverse(y):
    return np.digitize(y, x_bins)


def forward(x):
    y = x / (len(x_bins) * 10)
    return y


def append_scatter_params(df: pd.DataFrame, params: dict) -> pd.DataFrame:
    dfparams = pd.DataFrame.from_dict(params)
    # repeat values for all variants
    dfparams = pd.concat([dfparams] * df.shape[0], ignore_index=True)
    dfparams.set_index(df.index, inplace=True)
    # append plotting params to the points coordinates
    df = pd.concat([df, dfparams], axis=1)
    return df


def binned_countplot(df: pd.DataFrame, figpath: str):
    """
    :param df: columns needed:
        af_bin: dtype = float
        dataset: dtype = object
    :param figpath: path to directory where the figure is savec
    """
    with sns.axes_style("whitegrid"):
        g = sns.FacetGrid(data=df, col="dataset", margin_titles=True, height=figsize, aspect=1.5)
    g.map(sns.countplot, 'af_bin', color='0.3')
    # g.map(sns.distplot, 'af_info',  bins=x_bins, hist=True, kde=False, rug=False, color='0.3')  # not readable
    g.set_axis_labels("Estimated alternate allele frequency in main population", "Markers count")
    g.fig.subplots_adjust(wspace=.08, hspace=.02)
    g.despine(left=True)
    g.savefig(os.path.join(figpath,
                           'eda_pooled_gt_counts.pdf'))


def binned_striped_violin(df: pd.DataFrame, variable: str, ylabel: str, figpath: str):
    """
        :param df: columns needed:
            'af_bin': dtype = float
            'af_info': dtype = float
            'dataset': dtype = object
        :param variable: df's column to plot on y-axis, either 'missing_rate' or 'aaf'
        :param ylabel: label for the y-axis
        :param figpath: path to directory where the figure is savec
    """
    with sns.axes_style("whitegrid"):
        g = sns.FacetGrid(data=df, col="dataset", margin_titles=True, height=figsize, aspect=1.5)
    g.map(sns.stripplot, 'af_bin', variable,
          color=sns.color_palette("muted")[0], size=2.5, edgecolor='w', linewidth=0.2)
    g.map(sns.violinplot, 'af_bin', variable,
          bw='silverman', cut=1, width=1.0, inner=None, color=".8", linewidth=0.5)
    g.set_axis_labels("Estimated alternate allele frequency in main population", ylabel)
    if variable == 'aaf':
        g.set(yticks=lab_bins)
        g.set_yticklabels(spaced2_labs, fontsize=labfont)
    if variable == 'missing_rate':
        g.set(yticks=np.arange(0.0, 1.0, 0.1).round(decimals=1))
        g.set_yticklabels(np.arange(0.0, 1.0, 0.1).round(decimals=1), fontsize=labfont)
    g.fig.subplots_adjust(wspace=.08, hspace=.02)
    g.savefig(os.path.join(figpath,
                           'eda_pooled_gt_violin_{}.pdf'.format(variable)))


def binned_striped_box(df: pd.DataFrame, variable: str, ylabel: str, figpath: str):
    """
        :param df: columns needed:
            'af_bin': dtype = float
            'af_info': dtype = float
            'dataset': dtype = object
        :param variable: df's column to plot on y-axis, either 'missing_rate' or 'aaf'
        :param ylabel: label for the y-axis
        :param figpath: path to directory where the figure is savec
    """
    with sns.axes_style("whitegrid"):
        g = sns.FacetGrid(data=df, col="dataset", margin_titles=True, height=figsize, aspect=1.5)
    g.map(sns.stripplot, 'af_bin', variable,
          color=sns.color_palette("muted")[0], size=2.5, edgecolor='w', linewidth=0.2)
    g.map(sns.boxplot, 'af_bin', variable,
          whis=np.inf, color=".8", linewidth=0.5, fliersize=0.0)
    g.set_axis_labels("Estimated alternate allele frequency in main population", ylabel)
    if variable == 'aaf':
        g.set(yticks=lab_bins)
        g.set_yticklabels(spaced2_labs, fontsize=labfont)
    if variable == 'missing_rate':
        g.set(yticks=np.arange(0.0, 1.0, 0.1).round(decimals=1))
        g.set_yticklabels(np.arange(0.0, 1.0, 0.1).round(decimals=1), fontsize=labfont)
    g.fig.subplots_adjust(wspace=.08, hspace=.02)
    g.savefig(os.path.join(figpath,
                           'eda_pooled_gt_box_{}.pdf'.format(variable)))


def scatter_joint(df: pd.DataFrame, variable: str, ylabel: str, figtitle: str, figpath: str):
    """
        :param df: columns needed:
            'af_info': dtype = float
        :param variable: df's column to plot on y-axis
        :param ylabel: label for the y-axis
        :param figtitle: suptitle for the JointGrid
        :param figpath: path to directory where the figure is savec
    """
    sns.set(style="white", color_codes=True)
    g = sns.jointplot("af_info", variable, data=df,  # figure data
                      kind="scatter", color='0.5', s=10, edgecolor="w", linewidth=1,  # scatterplot params
                      marginal_kws=dict(bins=x_bins, hist=True, color='0.3'),  # distplot params
                      height=figsize, space=0, ratio=3)  # figure layout params
    g.set_axis_labels("Estimated alternate allele frequency in main population",
                      ylabel)
    plt.suptitle('data = {}'.format(figtitle))
    g.fig.subplots_adjust(wspace=.08, hspace=.02)
    g.savefig(os.path.join(figpath,
                           'eda_{}_gt_joint_{}.pdf'.format(figtitle, variable)))


# Drawing

df1['af_bin'] = pd.cut(df1.af_info, x_bins, labels=lab_bins)
df2['af_bin'] = pd.cut(df2.af_info, x_bins, labels=lab_bins)

# Plots 0: Missing rate in data
df10 = append_scatter_params(df1[['af_info', 'af_bin', 'missing_rate']], {'dataset': ['non-pooled']})
df20 = append_scatter_params(df2[['af_info', 'af_bin', 'missing_rate']], {'dataset': ['pooled']})
dfmiss = pd.concat([df10, df20], axis=0)
dfbin = dfmiss[['af_bin', 'dataset']].reset_index()
dfbin.astype(object, copy=False)

# countplot
binned_countplot(dfbin,
                 os.path.dirname(paths['gt']['pooled']))

# striped violin
binned_striped_violin(dfmiss,
                      'missing_rate',
                      "Genotype missing rate",
                      os.path.dirname(paths['gt']['pooled']))

# striped box
binned_striped_box(dfmiss,
                   'missing_rate',
                   "Genotype missing rate",
                   os.path.dirname(paths['gt']['pooled']))

# with sns.axes_style("whitegrid"):
#     g = sns.FacetGrid(data=dfmiss, col="dataset", margin_titles=True, height=4.5)
# g.map(sns.stripplot, 'af_bin', 'missing_rate',
#       color=sns.color_palette("muted")[0], size=2.5, edgecolor='w', linewidth=0.2)
# g.map(sns.boxenplot, 'af_bin', 'missing_rate', color='0.8', linewidth=0.2)
# g.set_axis_labels("Estimated alternate allele frequency in main population", "Genotype missing rate")
# # g.set(yticks=lab_bins) # yticks=np.arange(0.0, 1.0, 0.1))
# g.fig.subplots_adjust(wspace=.08, hspace=.02)
# plt.savefig(os.path.join(os.path.dirname(paths['gt']['pooled']), 'eda_pooled_gt_missing_boxenplot.pdf'))


# Plots 1: AAF in data
df11 = append_scatter_params(df1[['af_info', 'af_bin', 'aaf']], {'dataset': ['non-pooled']})
df21 = append_scatter_params(df2[['af_info', 'af_bin', 'aaf']], {'dataset': ['pooled']})
dfaaf = pd.concat([df11, df21], axis=0)

# jointplot
scatter_joint(df11,
              'aaf',
              "Computed alternate allele frequency in the study population",
              'non-pooled',
              os.path.dirname(paths['gt']['pooled']))

scatter_joint(df21,
              'aaf',
              "Computed alternate allele frequency in the study population",
              'pooled',
              os.path.dirname(paths['gt']['pooled']))

# striped violin
binned_striped_violin(dfaaf,
                      'aaf',
                      "Computed alternate allele frequency in the study population",
                      os.path.dirname(paths['gt']['pooled']))

# striped box
binned_striped_box(dfaaf,
                   'aaf',
                   "Computed alternate allele frequency in the study population",
                   os.path.dirname(paths['gt']['pooled']))


boxw = 0.75
plt.rcParams["figure.figsize"] = [figsize*2, figsize*boxw*1]
fig, axes = plt.subplots(1, 2)  # , sharey=True)

# Plot 20: genotypes proportions and HWE in true data
params200 = {'data': ['heterozygous']}
df200 = append_scatter_params(df3[['af_info', 'het_rate']], params200)
params201 = {'data': ['homozygous_alt']}
df201 = append_scatter_params(df3[['af_info', 'hom_alt_rate']], params201)
params202 = {'data': ['homozygous_ref']}
df202 = append_scatter_params(df3[['af_info', 'hom_ref_rate']], params202)

df20 = pd.concat([df200, df201, df202], axis=0)
df20 = pd.melt(df20, value_vars=['het_rate', 'hom_alt_rate', 'hom_ref_rate'], var_name='data', ignore_index=False)
df20 = df20.join(df3.af_info, how='inner')
df20.dropna(inplace=True)
df20.sort_values(by=['af_info', 'data'], inplace=True)
sns.scatterplot(data=df20, x='af_info', y='value', hue='data', s=pts_sz, ax=axes[0], linewidth=0.1, palette="GnBu_d")

hwe20 = np.linspace(0.0, 1.0, num=100)
rr20 = (1 - hwe20)**2
ra20 = 2 * (1 - hwe20) * hwe20
aa20 = hwe20**2
sns.lineplot(x=hwe20, y=aa20, color=sns.color_palette("GnBu_d")[2], ax=axes[0], label='hom_alt under HWE')
sns.lineplot(x=hwe20, y=ra20, color=sns.color_palette("GnBu_d")[0], ax=axes[0], label='het under HWE')
sns.lineplot(x=hwe20, y=rr20, color=sns.color_palette("GnBu_d")[4], ax=axes[0], label='hom_ref under HWE')

axes[0].set_xlabel("Estimated alternate allele frequency in main population")
axes[0].set_ylim(-0.01, 1.05)
axes[0].set_ylabel('Genotype proportion')
axes[0].set_title('data = non-pooled')
box = axes[0].get_position()
axes[0].set_position([box.x0, box.y0, box.width * boxw, box.height])  # resize position
axes[0].get_legend().set_visible(False)  # disable legend


# Plot 21: genotypes proportions and HWE in pooled data
params210 = {'data': ['heterozygous']}
df210 = append_scatter_params(df4[['af_info', 'het_rate']], params210)
params211 = {'data': ['homozygous_alt']}
df211 = append_scatter_params(df4[['af_info', 'hom_alt_rate']], params211)
params212 = {'data': ['homozygous_ref']}
df212 = append_scatter_params(df4[['af_info', 'hom_ref_rate']], params212)

df21 = pd.concat([df210, df211, df212], axis=0)
df21 = pd.melt(df21, value_vars=['het_rate', 'hom_alt_rate', 'hom_ref_rate'], var_name='data', ignore_index=False)
df21 = df21.join(df4.af_info, how='inner')
df21.dropna(inplace=True)
df21.sort_values(by=['af_info', 'data'], inplace=True)
sns.scatterplot(data=df21, x='af_info', y='value', hue='data', s=pts_sz, ax=axes[1], linewidth=0.1, palette="GnBu_d")

hwe10 = np.linspace(0.0, 1.0, num=100)
rr10 = (1 - hwe10)**2
ra10 = 2 * (1 - hwe10) * hwe10
aa10 = hwe10**2
sns.lineplot(x=hwe10, y=aa10, color=sns.color_palette("GnBu_d")[2], ax=axes[1], label='hom_alt under HWE')
sns.lineplot(x=hwe10, y=ra10, color=sns.color_palette("GnBu_d")[0], ax=axes[1], label='het under HWE')
sns.lineplot(x=hwe10, y=rr10, color=sns.color_palette("GnBu_d")[4], ax=axes[1], label='hom_ref under HWE')

axes[1].set_xlabel("Estimated alternate allele frequency in main population")
axes[1].set_ylim(-0.01, 1.05)
axes[1].set_ylabel('Genotype proportion')
axes[1].set_title('data = pooled')
box = axes[1].get_position()
axes[1].set_position([box.x0 * 0.85, box.y0, box.width * boxw, box.height])  # resize position
axes[1].legend(loc='center right', bbox_to_anchor=(1 + boxw, 0.5), ncol=1)  # Put a legend to the right side

plt.savefig(os.path.join(os.path.dirname(paths['gt']['pooled']), 'eda_genos_hwe.pdf'))
