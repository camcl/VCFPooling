import pandas as pd

from src.VCFPooling.poolSNPs import parameters as prm
from src.VCFPooling.poolSNPs import utils
from persotools.files import *
from persotools.debugging import *

dbg = MyPrintClass(True)

"""
Compute the imputation errors.
Write to csv
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


