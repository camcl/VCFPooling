from cyvcf2 import VCF
from scipy.stats import *
import numpy as np
import pandas as pd

from scripts.VCFPooling.poolSNPs.alleles import alleles_tools as alltls
from scripts.VCFPooling.poolSNPs.alleles import alleles_plots as allplt
from scripts.VCFPooling.poolSNPs import parameters as prm
from persotools.files import *
from persotools.debugging import *

dbg = MyPrintClass(True)

"""
Utils for data sets processing
NumPy or pandas DataFrame data sets
"""


def pgcd(a, b):
    """Computes the biggest common divider"""
    while b != 0:
        r = a % b
        a, b = b, r
    return a


def ppcm(a, b):
    """Computes the smallest common multiple"""
    if (a == 0) or (b == 0):
        return 0
    else:
        return (a*b)//pgcd(a, b)


def normalize(v):
    norm = np.linalg.norm(v)
    if norm == 0:
        return v
    else:
        return v / norm


#TODO: PEP8 refactoring
def sort_datasets(dflist: list, groups: list, df_aaf: pd.DataFrame) -> list:
    out = []
    df_aaf.reset_index(drop=True, inplace=True)
    for dfset in dflist:
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


