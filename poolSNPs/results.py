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
from scripts.poolSNPs import parameters as prm
from persotools.files import *


"""
Compute the imputation errors.
Vertically: read chunks of 1000 markers.
Plot different analyses.
"""

### TOOLS


def get_pop():
    df_pop = pd.read_csv('20130606_sample_info.csv',
                      sep=',',
                      header=0,
                      usecols=[0, 1])
    df_pop.sort_values(['Sample', 'Population'], axis=0, inplace=True)
    return df_pop


def sort_datasets(args, groups, df_maf):
    out = []
    for dfset in args:
        # Sort samples by population
        dfset.sort_index(axis=1, inplace=True)
        dfset.columns = groups
        dfset.sort_index(level='Population', axis=1, inplace=True)
        # Sort variants by MAF
        dfset.sort_index(axis=0, inplace=True) # sort by id
        dfset.reset_index(drop=True, inplace=True)
        dfset['maf'] = df_maf['maf']
        dfset['maf_inter'] = df_maf['maf_inter']
        dfset.set_index([df_maf['id'], 'maf', 'maf_inter'], drop=True, append=False, inplace=True)  # replace idx with multiidx (id sorted)
        dfset.sort_index(level=['maf_inter', 'maf'], axis=0, inplace=True)
        dfset.reset_index(drop=True, inplace=True)
        out.append(dfset)
    return out


def make_index(raw_data):
    src = prm.SOURCE
    group = get_pop()
    df_maf = alleles_tools.get_maf(src)
    samples = pd.DataFrame(VCF(raw_data).samples, columns=['Sample'])
    df_pop = group.merge(samples, on='Sample', how='inner')

    maf_idx = pd.MultiIndex.from_arrays(list(zip(*df_maf.values)), names=['id', 'maf', 'maf_inter'])
    pop_idx = pd.MultiIndex.from_arrays(list(zip(*df_pop.values)), names=['Sample', 'Population'])

    return maf_idx, pop_idx


if __name__ == '__main__':
    ### DATA
    print('Configure parameters'.ljust(80, '.'))
    sorting = True # sort data sets by MAF and population values
    chk_sz = prm.CHK_SZ
    os.chdir(prm.WD)

    ALL = prm.SOURCE
    RAW = prm.RAW['b1i'] + '.vcf.gz'
    MISS = 'IMP.chr20.missing.beagle2.corr.vcf.gz'
    POOL = 'IMP.chr20.pooled.beagle2.corr.vcf.gz'
    #MISSPOOL = VCF('IMP.chr20.missing.pooled.beagle2.corr.vcf.gz')

    print('Raw samples and index labels'.ljust(80, '.'))
    r_size, ms_size = len(VCF(RAW).samples), len(VCF(MISS).samples)

    # Load MAFs
    mafs = alleles_tools.get_maf('ALL.chr20.snps.gt.chunk{}.vcf.gz'.format(str(chk_sz)))

    # Create line-index and column-index (multi-indexing)
    maf_idx, pop_idx = make_index(RAW)

    ### PROCESSING
    print('Build datasets'.ljust(80, '.'))
    raw0, raw1 = alleles_tools.vcf2dframe(VCF(RAW), r_size)
    miss0, miss1 = alleles_tools.vcf2dframe(VCF(MISS), ms_size)
    # misspool0, misspool1 = alleles_tools.vcf2dframe(MISSPOOL, ms_size)
    pool0, pool1 = alleles_tools.vcf2dframe(VCF(POOL), ms_size)
    datasets = [raw0, raw1, miss0, miss1, pool0, pool1]

    if sorting:
        raw0, raw1, miss0, miss1, pool0, pool1 = sort_datasets(datasets, pop_idx, mafs)

    samples = VCF(RAW).samples
    variants = raw0.index.tolist()

    # NumPy juxtapos for error computation
    raw = np.stack([raw0.values, raw1.values], axis=-1)
    miss = np.stack([miss0.values, miss1.values], axis=-1)
    pool = np.stack([pool0.values, pool1.values], axis=-1)

    set_errors = alleles_tools.compute_imp_err(['pooled', 'missing'],
                                 [pool, miss],
                                 raw,
                                 maf_idx,
                                 pop_idx,
                                 sorting)

    ### PLOT BOXPLOTS AND DENSITY CURVES OF IMPUTATION ERRORS
    # alleles_plots.boxplot_densities(set_errors)

    ### FOCUS ON THE LOWEST MAF MARKERS
    import sys
    #sys.stdout = open('.log', 'w')

    alleles_plots.multiplot_err_het(set_errors)
    
    ### ERROR VS. HETEROZIGOSITY/MISSING ALLELES
    for k, err in set_errors.items():
        reffile = prm.RAW['imp']
        alleles_plots.plot_err_vs_het(k,
                                      err,
                                      reffile,
                                      low_maf=False,
                                      save=True)
        alleles_plots.plot_err_vs_miss(k,
                                       err)

    ### MAF BEHAVIOR
    alleles_plots.plot_mafs(set_errors, ALL, boxplot=False, lineplot=True)

    ### HEATMAPS OF IMPUTATION ERROR PER SITE PER SAMPLE
    # for k, err in set_errors.items():
    #     pass
    #     print('Plotting {} errors...'.format(k))
    #     alleles_plots.plot_heat_map(err,
    #                   k.upper(),
    #                   [16, 4],
    #                   sorting,
    #                   title="Errors in imputation (Beagle 4.1): RAW data vs. {} data",
    #                   rightYaxis=False)
    #     alleles_plots.plot_heat_map(err.query('maf_inter >= 2'),
    #                   k.upper() + '.low_maf',
    #                   [8, 6],
    #                   sorting,
    #                   title="Errors in imputation (Beagle 4.1) of lowest-MAF variants: RAW data vs. {} data",
    #                   rightYaxis=False)

    # print('Correlation\'s coeff'.ljust(80, '.'))
    # for k, sp in spearmans.items():
    #     print(k, '\r')
    #     print(np.array(sp[0][0,:,:]), np.array(sp[1][0,:,:])) # rho mtx for each axis

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