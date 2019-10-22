from cyvcf2 import VCF
from scipy.stats import *
import numpy as np
import pandas as pd

from scripts.VCFPooling.poolSNPs.alleles import alleles_tools as alltls
from scripts.VCFPooling.poolSNPs.alleles import alleles_plots as allplt
from scripts.VCFPooling.poolSNPs import parameters as prm
from scripts.VCFPooling.poolSNPs import utils
from persotools.files import *
from persotools.debugging import *

dbg = MyPrintClass(True)

"""
Compute the imputation errors.
Vertically: read chunks of 1000 markers.
Plot different analyses.
"""

### TOOLS


def sort_datasets(*x):
    return utils.sort_datasets(*x)


if __name__ == '__main__':
    ### DATA
    print('Load parameters'.ljust(80, '.'))
    sorting = True  # sort data sets by AAF and population values
    subsample = prm.SUBSET
    chk_sz = prm.CHK_SZ

    cd = os.path.join(prm.WD, 'gt', 'stratified', 'all_snps_all_samples')
    os.chdir(cd)

    ALL = os.path.join(prm.WD, 'gt', 'stratified', prm.CHKFILE)
    RAW = os.path.join(prm.WD, 'gt', 'stratified',
                       prm.RAW['b1i'] + '.vcf.gz')
    MISS = os.path.join(prm.WD, 'gt', 'stratified',
                        'all_snps_all_samples', prm.MISSING['gtonly'] + '.vcf.gz')
    POOL = os.path.join(prm.WD, 'gt', 'stratified',
                        'all_snps_all_samples', prm.POOLED['gtonly'] + '.vcf.gz')

    vcfpathdic = {'all': ALL,
                  'raw': RAW,
                  'missing': MISS,
                  'pooled': POOL}

    print('Raw samples and index labels'.ljust(80, '.'))
    popsize = prm.NB_IMP

    set_errors = {}
    for k in ['pooled', 'missing']:
        dfmse = pd.read_csv(k + '{}.chunk{}.csv'.format('.sorted' if sorting else '', prm.CHK_SZ),
                            sep='\t',
                            encoding='utf-8',
                            index_col=0,
                            usecols=[0] + list(range(3, popsize + 3)),
                            skiprows=[0, 1]
                            )
        set_errors[k] = dfmse

    '''

    ### PLOT BOXPLOTS AND DENSITY CURVES OF IMPUTATION ERRORS
    # allplt.boxplot_densities(set_errors)

    # allplt.multiplot_err_het(set_errors)

    ### ERROR VS. HETEROZIGOSITY/MISSING ALLELES
    # for k, err in set_errors.items():
    #     reffile = prm.RAW['imp']
    #     allplt.plot_err_vs_het(k,
    #                                   err,
    #                                   reffile,
    #                                   low_aaf=False,
    #                                   save=True)
    #     allplt.plot_err_vs_miss(k,
    #                                    err)
    '''

    ### AAF BEHAVIOR
    # allplt.plot_aaf_evol(set_errors, ALL)
    for setkey, setval in set_errors.items():
        allplt.plot_aaf_twist(setkey, setval, vcfpathdic[setkey])

    # TODO: correlation between aaf and aaf imputed. Move in a alleles.metrics file
