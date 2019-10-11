from cyvcf2 import VCF
from scipy.stats import *
import numpy as np
import pandas as pd

from scripts.poolSNPs.alleles import alleles_tools as alltls
from scripts.poolSNPs.alleles import alleles_plots as allplt
from scripts.poolSNPs import parameters as prm
from persotools.files import *
from persotools.debugging import *

dbg = MyPrintClass(True)

"""
Utils for data sets processing
NumPy or pandas DataFramde data sets
"""

### TOOLS


def sort_datasets(args, groups, df_aaf):
    out = []
    df_aaf.reset_index(drop=True, inplace=True)
    for dfset in args:
        # Sort samples by population
        dfset.sort_index(axis=1, inplace=True)
        dfset.columns = groups
        dfset.sort_index(level='Population', axis=1, inplace=True)
        # Sort variants by AAF
        dfset.sort_index(axis=0, inplace=True) # sort by id
        dfset.reset_index(drop=True, inplace=True)
        dfset['af_info'] = df_aaf['af_info']
        dfset['aaf_bin'] = df_aaf['aaf_bin']
        dfset.set_index([df_aaf['id'], 'af_info', 'aaf_bin'], drop=True, append=False, inplace=True)  # replace idx with multiidx (id sorted)
        dfset.sort_index(level=['aaf_bin', 'af_info'], axis=0, inplace=True)
        dfset.reset_index(drop=True, inplace=True)
        out.append(dfset)
    return out