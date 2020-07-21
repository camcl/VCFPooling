from scripts.VCFPooling.poolSNPs import BCFstatsUtils as bcfstats
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

"""
This code reproduces part of the customized metrics in quality.py.
Data used is generated via the `bcftools stats` command and stored in a text-like file.
The metrics are suitable only for VCF files having a GT field.
"""

def get_counts():
    """
    Matches and mismatches counts from bcfstats data frame
    """
    pass


def get_overall_nrd(filestats: str) -> pd.DataFrame:
    """
    :param filestats: path to file
    # NRD and discordance is calculated as follows:
    #   m .. number of matches
    #   x .. number of mismatches
    #   NRD = (xRR + xRA + xAA) / (xRR + xRA + xAA + mRA + mAA)
    #   RR discordance = xRR / (xRR + mRR)
    #   RA discordance = xRA / (xRA + mRA)
    #   AA discordance = xAA / (xAA + mAA)
    """
    nrds = bcfstats.get_table_dataframe(filestats, "NRDs")

    return nrds


def get_bin_nrd(filestats: str) -> pd.DataFrame:
    """
    Calculates NRD for each bin of AF
    :param filestats: path to file
    # NRD and discordance is calculated as follows:
    #   m .. number of matches
    #   x .. number of mismatches
    #   NRD = (xRR + xRA + xAA) / (xRR + xRA + xAA + mRA + mAA)
    #   RR discordance = xRR / (xRR + mRR)
    #   RA discordance = xRA / (xRA + mRA)
    #   AA discordance = xAA / (xAA + mAA)
    """
    genos = {0: 'RR Hom', 1: 'RA Het', 2: 'AA Hom'}
    gcsaf = bcfstats.get_table_dataframe(filestats, "GCsAF")
    gcsaf = gcsaf.astype(float)
    gcsaf['tot matches'] = gcsaf['RR Hom matches'] + gcsaf['RA Het matches'] + gcsaf['AA Hom matches']
    gcsaf['tot mismatches'] = gcsaf['RR Hom mismatches'] + gcsaf['RA Het mismatches'] + gcsaf['AA Hom mismatches']

    gcsaf['NRD'] = (gcsaf['RR Hom mismatches'] + gcsaf['RA Het mismatches'] + gcsaf['AA Hom mismatches']) /\
                   (gcsaf['RR Hom mismatches'] + gcsaf['RA Het mismatches'] + gcsaf['AA Hom mismatches'] + gcsaf['RA Het matches'] + gcsaf['AA Hom matches'])
    for geno in genos.values():
        gcsaf['{} discordance'.format(geno)] = gcsaf['{} mismatches'.format(geno)] / (gcsaf['{} mismatches'.format(geno)] + gcsaf['{} matches'.format(geno)])

    return gcsaf


def get_bin_maf_nrd(filestats: str) -> pd.DataFrame:
    """
    Calculates NRD for each bin of AF
    :param filestats: path to file
    # NRD and discordance is calculated as follows:
    #   m .. number of matches
    #   x .. number of mismatches
    #   NRD = (xRR + xRA + xAA) / (xRR + xRA + xAA + mRA + mAA)
    #   RR discordance = xRR / (xRR + mRR)
    #   RA discordance = xRA / (xRA + mRA)
    #   AA discordance = xAA / (xAA + mAA)
    """
    genos = {0: 'RR Hom', 1: 'RA Het', 2: 'AA Hom'}
    gcsaf = bcfstats.get_table_dataframe(filestats, "GCsAF")
    gcsaf = gcsaf.astype(float)

    # Create MAF values from AF and aggregate the counts for MAF
    mafconverter = lambda x: x if x <= 0.5 else 1 - x
    gcsaf['minor allele frequency'] = gcsaf['allele frequency'].apply(mafconverter).round(decimals=2)
    gcsaf.drop(labels=['allele frequency', 'dosage r-squared'], axis=1, inplace=True)
    gpby = gcsaf.groupby(['minor allele frequency'], axis=0)
    gcsaf = gpby.sum()

    # Compute number of (mis)matches
    gcsaf['tot matches'] = gcsaf['RR Hom matches'] + gcsaf['RA Het matches'] + gcsaf['AA Hom matches']
    gcsaf['tot mismatches'] = gcsaf['RR Hom mismatches'] + gcsaf['RA Het mismatches'] + gcsaf['AA Hom mismatches']

    gcsaf['NRD'] = (gcsaf['RR Hom mismatches'] + gcsaf['RA Het mismatches'] + gcsaf['AA Hom mismatches']) / \
                   (gcsaf['RR Hom mismatches'] + gcsaf['RA Het mismatches'] + gcsaf['AA Hom mismatches'] + gcsaf[
                       'RA Het matches'] + gcsaf['AA Hom matches'])
    for geno in genos.values():
        gcsaf['{} discordance'.format(geno)] = gcsaf['{} mismatches'.format(geno)] / (
                    gcsaf['{} mismatches'.format(geno)] + gcsaf['{} matches'.format(geno)])

    return gcsaf


def get_bin_maf_concordance(filestats: str) -> pd.DataFrame:
    """
    Calculates NRD for each bin of AF
    :param filestats: path to file
    # Concordance is calculated as follows:
    #   m .. number of matches
    #   x .. number of mismatches
    #   concordance GC = (mxRR + mRA + mAA) / (xRR + xRA + xAA + mRA + mAA + mRR)
    #   RR concordance = mRR / (xRR + mRR)
    #   RA concordance = mRA / (xRA + mRA)
    #   AA concordance = mAA / (xAA + mAA)
    """
    genos = {0: 'RR Hom', 1: 'RA Het', 2: 'AA Hom'}
    gcsaf = bcfstats.get_table_dataframe(filestats, "GCsAF")
    gcsaf = gcsaf.astype(float)

    # Create MAF values from AF and aggregate the counts for MAF
    mafconverter = lambda x: x if x <= 0.5 else 1 - x
    gcsaf['minor allele frequency'] = gcsaf['allele frequency'].apply(mafconverter).round(decimals=2)
    gcsaf.drop(labels=['allele frequency', 'dosage r-squared'], axis=1, inplace=True)
    gpby = gcsaf.groupby(['minor allele frequency'], axis=0)
    gcsaf = gpby.sum()

    # Compute number of (mis)matches
    gcsaf['tot matches'] = gcsaf['RR Hom matches'] + gcsaf['RA Het matches'] + gcsaf['AA Hom matches']
    gcsaf['tot mismatches'] = gcsaf['RR Hom mismatches'] + gcsaf['RA Het mismatches'] + gcsaf['AA Hom mismatches']

    gcsaf['GC'] = gcsaf['tot matches'] / (gcsaf['tot mismatches'] + gcsaf['tot matches'])
    for geno in genos.values():
        gcsaf['{} concordance'.format(geno)] = gcsaf['{} matches'.format(geno)] / (
                    gcsaf['{} mismatches'.format(geno)] + gcsaf['{} matches'.format(geno)])

    return gcsaf


def get_bin_accuracy(filestats: str) -> pd.DataFrame:
    """
    :param filestats: path to file
    ratio tp / (tp + fp + tn + fn) for each class
    """
    genos = {0: 'RR Hom', 1: 'RA Het', 2: 'AA Hom'}
    gcsaf = bcfstats.get_table_dataframe(filestats, "GCsAF")
    gcsaf = gcsaf.astype(float)

    for geno in genos.values():
        gcsaf['{} accuracy'.format(geno)] = gcsaf['{} matches'.format(geno)] / gcsaf['number of genotypes']

    return gcsaf

def get_bin_precision(filestats: str) -> pd.DataFrame:
    """
    :param filestats: path to file
    ratio tp / (tp + fp) where tp is the number of true positives and fp the number of false positives
    """
    # this cannot be calculated: I do not know which one are the false positives?
    genos = {0: 'RR Hom', 1: 'RA Het', 2: 'AA Hom'}
    gcsaf = bcfstats.get_table_dataframe(filestats, "GCsAF")
    gcsaf = gcsaf.astype(float)

    for geno in genos.values():
        pass
        #gcsaf['{} accuracy'.format(geno)] = gcsaf['{} matches'.format(geno)] / gcsaf['number of genotypes']
    # or?
    # gcsaf['RR Hom precision'] = gcsaf['RR Hom matches'] /
    # (gcsaf['RA Het mismatches'] + gcsaf['AA Hom mismatches'] - gcsaf['RR Hom mismatches'] + gcsaf['RR Hom matches'])

    return None


def get_bin_recall(filestats: str) -> pd.DataFrame:
    """
    :param filestats: path to file
    ratio tp / (tp + fn) for each class
    """
    # this is the per-class concordance
    genos = {0: 'RR Hom', 1: 'RA Het', 2: 'AA Hom'}
    gcsaf = bcfstats.get_table_dataframe(filestats, "GCsAF")
    gcsaf = gcsaf.astype(float)

    for geno in genos.values():
        gcsaf['{} recall'.format(geno)] = gcsaf['{} matches'.format(geno)] / (gcsaf['{} mismatches'.format(geno)] + gcsaf['{} matches'.format(geno)])

    return gcsaf