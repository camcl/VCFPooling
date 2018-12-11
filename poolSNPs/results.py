from cyvcf2 import VCF
import os
import numpy as np
import pandas as pd
import subprocess
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.stats import *
import math
from scripts.poolSNPs.alleles import alleles_tools
from scripts.poolSNPs.alleles import alleles_plots


"""
Compute the imputation errors.
Vertically: read chunks of 1000 markers.
Plot different analysis.
"""
# TODO: Unphased poolning och prephasa data: kör Beagle med eller uten referens data set?
# TODO: Köra ny experiment med 100 000 markörer
# TODO: try to characterize behav# TODO: try to characterize behavior of markers

### TOOLS


def get_pop():
    df_pop = pd.read_csv('20130606_sample_info.csv',
                      sep=',',
                      header=0,
                      usecols=[0, 1])
    df_pop.sort_values(['Sample', 'Population'], axis=0, inplace=True)
    return df_pop


def sort_datasets(args, pop_idx, MAF):
    out = []
    for dfset in args:
        # Sort samples by population
        dfset.sort_index(axis=1, inplace=True)
        dfset.columns = pop_idx
        dfset.sort_index(level='Population', axis=1, inplace=True)
        # Sort variants by MAF
        dfset.sort_index(axis=0, inplace=True) # sort by id
        dfset.reset_index(drop=True, inplace=True)
        dfset['maf'] = MAF['maf']
        dfset['maf_inter'] = MAF['maf_inter']
        dfset.set_index([MAF['id'], 'maf', 'maf_inter'], drop=True, append=False, inplace=True)  # replace idx with multiidx (id sorted)
        dfset.sort_index(level=['maf_inter', 'maf'], axis=0, inplace=True)
        dfset.reset_index(drop=True, inplace=True)
        out.append(dfset)
    return out


if __name__ == '__main__':
    ### DATA
    print('Configure parameters'.ljust(80, '.'))
    sorting = True # sort data sets by MAF and population values
    os.chdir('/home/camilleclouard/PycharmProjects/1000Genomes/data/tests-beagle')

    RAW = VCF('IMP.chr20.beagle1.vcf.gz')
    MISS = VCF('IMP.chr20.missing.beagle2.corr.vcf.gz')
    POOL = VCF('IMP.chr20.pooled.beagle2.corr.vcf.gz')
    MISSPOOL = VCF('IMP.chr20.missing.pooled.beagle2.corr.vcf.gz')

    print('Raw samples and index labels'.ljust(80, '.'))
    r_size, ms_size = len(RAW.samples), len(MISS.samples)
    #TODO: implement 1000-sized chunking

    MAF = alleles_tools.get_maf('ALL.chr20.snps.gt.chunk.vcf.gz')
    maf_idx = pd.MultiIndex.from_arrays(list(zip(*MAF.values)), names=['id', 'maf', 'maf_inter'])
    POP = get_pop()
    samples = pd.DataFrame(RAW.samples, columns=['Sample'])
    pop = POP.merge(samples, on='Sample', how='inner')
    pop_idx = pd.MultiIndex.from_arrays(list(zip(*pop.values)), names=['Sample', 'Population'])

    ### PROCESSING
    print('Build datasets'.ljust(80, '.'))
    raw0, raw1 = alleles_tools.vcf2dframe(RAW, r_size)
    miss0, miss1 = alleles_tools.vcf2dframe(MISS, ms_size)
    misspool0, misspool1 = alleles_tools.vcf2dframe(MISSPOOL, ms_size)
    pool0, pool1 = alleles_tools.vcf2dframe(POOL, ms_size)
    datasets = [raw0, raw1, miss0, miss1, pool0, pool1, misspool0, misspool1]

    if sorting:
        raw0, raw1, miss0, miss1, pool0, pool1, misspool0, misspool1 = sort_datasets(datasets, pop_idx, MAF)

    samples = samples.values
    variants = raw0.index.tolist()

    raw = np.stack([raw0.values, raw1.values], axis=-1)
    miss = np.stack([miss0.values, miss1.values], axis=-1)
    misspool = np.stack([misspool0.values, misspool1.values], axis=-1)
    pool = np.stack([pool0.values, pool1.values], axis=-1)

    set_errors = alleles_tools.compute_imp_err(['pooled', 'missing', 'missing.pooled'],
                                 [pool, miss, misspool],
                                 raw,
                                 maf_idx,
                                 pop_idx,
                                 sorting)

    '''
    ### PLOT BOXPLOTS AND DENSITY CURVES OF IMPUTATION ERRORS
    pd.set_option('precision', 25)
    plt.rcParams["figure.figsize"] = [20, 20]
    plt.rcParams["figure.autolayout"] =  True
    errbox, axs = plt.subplots(3,2)

    pos = 0
    for k in ['pooled', 'missing', 'missing.pooled']:
        df = set_errors[k]['grid']
        maf_err = df.groupby(axis=0, level='maf_inter').mean()
        pop_err = df.groupby(axis=1, level='Population').mean()
        lev_err = df.applymap(lambda x: np.power(x,2)).apply(np.mean,
                                               axis=1).rename('log1p_mse').to_frame().applymap(math.log1p)
        lev_err.reset_index(drop=False, inplace=True)
        lev_err.drop(['id', 'maf'], axis=1, inplace=True)
        lev_err.boxplot(by='maf_inter',
                        ax=axs[pos, 0],
                        showmeans=True,
                        showfliers=False)
        plt.suptitle("")
        axs[pos, 0].set_title('log1p_mse for the {} data set'.format(k))
    
        cross_err = np.power(maf_err, 2).groupby(axis=1,
                                             level='Population').mean().applymap(math.log1p)
        print('\r\nNumber of missing genotypes in the "{}" dataset:'.format(k))
        nbmiss = subprocess.call("bcftools view -H IMP.chr20.{}.cfgt.vcf.gz | grep -o '\.\|\.' | wc -l".format(k),
                                 shell=True)
        print('Mean Error: {}'.format(set_errors[k]['mean']))
        # Plotting grouped errors
        # print('Log10-Errors per MAF\'s and population\'s category:\r', cross_err)
    
        axs[pos, 1].set_xlim(left=-1.0, right=1.0)
        pop_err.plot.density(ax=axs[pos, 1],
                             legend=False)
        # pop_err.loc[pop_err.index.levels[1].map(lambda x: x <= 0.05)].plot.density(ax=axs[pos, 1],
        #                      legend=False)
        pop_skw = skew(pop_err.values, axis=0, bias=True)
        axs[pos, 1].set_title('Error density distributions for the {} data set'.format(k))
        pos += 1
    handles, labels = axs[1,1].get_legend_handles_labels()
    detailed_labels = [pop + ', skewness = ' + str(sk) for pop, sk in zip(labels, pop_skw)]
    axs[1,1].legend(handles,
                    detailed_labels,
                    loc='center left',
                    bbox_to_anchor=(1.05, 0.5),
                    borderaxespad=0,
                    frameon=True)
    #plt.tight_layout()
    plt.savefig('log1p.mean.square.error.maf.level.png', dpi='figure')
    #plt.show()
    '''

    '''
    ### FOCUS ON THE LOWEST MAF MARKERS
    #TODO: boxplot for gps 2 & 3 + label stats values
    plt.rcParams["figure.figsize"] = [12, 6]
    low_maf, axm = plt.subplots(1,2)
    box = 0
    for k, dset in dict(zip(['missing', 'missing.pooled'], [miss, misspool])).items():
        df = set_errors[k]['grid']
        low_err = df.applymap(lambda x: np.power(x,2)).apply(np.mean,
                                               axis=1).rename('log1p_mse').to_frame().applymap(math.log1p)
        low_err.reset_index(drop=False, inplace=True)
        low_err.drop(['id', 'maf'], axis=1, inplace=True)
        low_err.set_index('maf_inter', inplace=True)
        low_err.query('maf_inter >= 2').boxplot(by='maf_inter',
                                                ax=axm[box],
                                                showmeans=True,
                                                showfliers=False)
        means = low_err.query('maf_inter >= 2').groupby(axis=0, level='maf_inter').mean()
        medians = low_err.query('maf_inter >= 2').median(axis=0, level='maf_inter')
        for mn, md in zip(means, medians):
            mpl.text.Annotation('mean = ' + str(mn), [mn.index, mn])
            mpl.text.Annotation('median = ' + str(md), [md.index, md])
        plt.suptitle("")
        axm[box].set_title('log1p_mse for the {} data set'.format(k))
        box += 1
    plt.savefig('log1p.mean.square.error.maf.low.png', dpi='figure')
    
    ### ERROR VS. HETEROZIGOSITY/MISSING ALLELES
    for k in set_errors.keys():
        break
        site_err = set_errors[k]['grid'].\
            applymap(lambda x: np.power(x,2)).\
            apply(np.mean, axis=1).\
            rename('log1p_mse').\
            to_frame().\
            applymap(math.log1p)
        site_err.reset_index(drop=False, inplace=True)
        site_err.drop(['maf_inter'], axis=1, inplace=True)
        site_err['het'] = site_err['id'].\
            map(per_site_heteroz('IMP.chr20.{}.cfgt.vcf.gz'.format(k))).\
            apply(lambda h: h*100/set_errors[k]['grid'].shape[1])
        site_err.plot.scatter('het', 'log1p_mse', c='maf', cmap=plt.cm.viridis)
        plt.axhline(y=set_errors[k]['log1p_mse'],
                    linestyle='dashed',
                    color='k',
                    label='mean error')
        plt.title('Imputation error vs. heterozygosity rate in {} data set'.format(k))
        # for s, x, y in zip(*[site_err[col].values for col in site_err.columns]):
        #     plt.text(x, y, s)
        #plt.show()
    '''

    ### MAF BEHAVIOR
    alleles_plots.plot_mafs(set_errors, boxplot=False, lineplot=True)

'''
    ### COMPUTE PEARSON'S CORRELATION COEFFICIENTS
    #TODO:
    # pearsons = {}
    # for k, dset in dict(zip(['miss', 'pool', 'misspool'], [miss, pool, misspool])).items():
    #     pearsons[k] = pearsonr(raw, dset)

    ### COMPUTE SPEARMAN'S CORRELATION COEFFICIENTS FOR VARIANTS AND SAMPLES
    spearmans = {}
    for k, err in set_errors.items():
        break
        spearmans[k] = [spearmanr(err, axis=0), spearmanr(err, axis=1)]
        #TODO: weird white lines on the plots?!
        alleles_plots.plot_heat_map(pd.DataFrame(spearmans[k][0].correlation),
                      k.upper() + '.spearman.lin',
                      [20, 20],
                      sorting,
                      title="Spearman's correlation between imputation's errors for {}",
                      legend="Correlation's coefficient")
        alleles_plots.plot_heat_map(pd.DataFrame(spearmans[k][1].correlation),
                      k.upper() + '.spearman.col',
                      [20, 20],
                      sorting,
                      title = "Spearman's correlation between imputation's errors for {}",
                      legend = "Correlation's coefficient")
'''
"""
    ### HEATMAPS OF IMPUTATION ERROR PER SITE PER SAMPLE
    for k, err in set_errors.items():
        print('Plotting {} errors...'.format(k))
        alleles_plots.plot_heat_map(err,
                      k.upper(),
                      [20, 10],
                      sorting,
                      title="Errors in imputation (Beagle 4.1): RAW data vs. {} data",
                      rightYaxis=True)
        alleles_plots.plot_heat_map(err.query('maf_inter >= 2'),
                      k.upper() + '.low_maf',
                      [20,3],
                      sorting,
                      title="Errors in imputation (Beagle 4.1) of lowest-MAF variants: RAW data vs. {} data",
                      rightYaxis=False)
"""

    # print('Correlation\'s coeff'.ljust(80, '.'))
    # for k, sp in spearmans.items():
    #     print(k, '\r')
    #     print(np.array(sp[0][0,:,:]), np.array(sp[1][0,:,:])) # rho mtx for each axis

