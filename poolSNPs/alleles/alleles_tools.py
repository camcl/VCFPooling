from cyvcf2 import VCF, Writer, Variant
import os
import numpy as np
import pandas as pd
import subprocess
import collections
from scipy.stats import *
import math
import scripts.poolSNPs.parameters as prm
from persotools.debugging import *
from persotools.files import *
import pickle
import inspect


"""
Compute the imputation errors.
Vertically: read chunks of 1000 markers.
Plot different analysis.
"""
# TODO: Unphased poolning och prephasa data: kör Beagle med eller uten referens data set?
# TODO: Köra ny experiment med 100 000 markörer
# TODO: try to characterize behav# TODO: try to characterize behavior of markers


### GENERAL TOOLS

def pgcd(a, b):
    """pgcd(a,b): calcul du 'Plus Grand Commun Diviseur' entre les 2 nombres entiers a et b"""
    while b != 0:
        r = a % b
        a, b = b, r
    return a


def ppcm(a, b):
    """ppcm(a,b): calcul du 'Plus Petit Commun Multiple' entre 2 nombres entiers a et b"""
    if (a == 0) or (b == 0):
        return 0
    else:
        return (a*b)//pgcd(a, b)


### SPECIFIC TOOLS

def per_site_sample_error(elem2x2):
    """
    Computes the imputation error as the distance between
    the true data set vs. imputed at one site for one sample.
    :param elem2x2: 2*2 array representing the pair of GTs to compare.
    Ex. 0|0  vs. 1|0 = [[0 0] [0 1]]
    :return: Z-norm. score for error: 0 = right, > 0 = wrong, max = 1
    """
    return np.array(abs(elem2x2[:2].sum() - elem2x2[2:].sum()))/2


def para_array(arr1, arr2):
    """
    Merge the two data sets to compare, adding 1 dimension.
    :param arr1: array of genotypes values
    :param arr2: array of genotypes values
    :return: array of genotypes values
    """
    return np.stack([arr1, arr2], axis=-1)


def per_element_error(arr4d, ax):
    """
    Computes imputation error row- or column-wise
    :param arr4d: data sets to compare
    :param ax: axis over which to apply the cumulative sum
    :return: 3d-array of error values from the compared data sets
    """
    arr = np.rollaxis(arr4d, ax)
    arr2d = np.reshape(arr, (arr.shape[0]*arr.shape[1], 2*2))
    errors = np.apply_along_axis(per_site_sample_error, 1, arr2d)
    return errors.reshape(arr.shape[0], arr.shape[1])


def per_axis_error(arr4d, ax):
    """
       Cumulates error row- or column-wise
       :param arr4d: data sets to compare
       :param ax: axis over which to apply the cumulative sum
       :return: 1d-array of error values along the ax-axis
       """
    return per_element_error(arr4d, ax).sum(axis=1)/arr4d.shape[ax]


def vcf2array(vcf_obj, size):
    """
    Throws the genotypes values of a VCF file into an array.
    :param vcf_obj: VCF file read with cyvcf2
    :return: 3d-array
    """
    vars = []
    arr = np.zeros((1, size, 2), dtype=int)
    for v in vcf_obj:
        vars.append(v.ID)
        arr = np.vstack((arr,
                         np.expand_dims(np.array(v.genotypes)[:, :-1],
                                        axis=0)))
    return arr[1:,:,:], vars


def vcf2dframe(vcf_obj, size, id='id'):
    """
    Throws the genotypes values of a VCF file into an array.
    :param vcf_obj: VCF file read with cyvcf2
    :return: two dataframes, one for each allele of the genotype
    """
    vars = []
    arr = np.zeros((1, size, 2), dtype=int)
    for v in vcf_obj:
        if id == 'id':
            vars.append(v.ID)
        if id == 'chrom:pos':
            vars.append(':'.join([str(v.CHROM), str(v.POS)]))
        arr = np.vstack((arr,
                         np.expand_dims(np.array(v.genotypes)[:, :-1],
                                        axis=0)))
    df0 = pd.DataFrame(arr[1:,:,0], index=vars, columns=vcf_obj.samples)
    df1 = pd.DataFrame(arr[1:, :, 1], index=vars, columns=vcf_obj.samples)
    return df0, df1


def convert_maf(x):
    try:
        x = float(x)
    except TypeError:
        x = np.nan
    finally:
        return x


def get_maf(vcf_raw, id='id'):
    # subprocess.run(['''bcftools query -f '%ID\t%AF\n' {0} > {1}/TMP.maf.csv'''.format(vcf_raw, os.getcwd())],
    #                 shell=True,
    #                 cwd='/home/camille/PycharmProjects/1000Genomes/data/tests-beagle')
    # df_maf = pd.read_csv('TMP.maf.csv',
    #                   sep='\t',
    #                   names=['id', 'maf'])
    dico = {'id': list(),
            'maf': list(), # alternate allele frequency estimated on the populations
            'aaf': list()} # alternate allele frequency estimated on the file
    if id == 'id':
        for var in VCF(vcf_raw):
            dico['id'].append(var.ID)
            dico['maf'].append(var.INFO['AF'])
            dico['aaf'].append(var.aaf)
    if id == 'chrom:pos':
        for var in VCF(vcf_raw):
            dico['id'].append(':'.join([str(var.CHROM), str(var.POS)]))
            dico['maf'].append(var.INFO['AF'])
            dico['aaf'].append(var.aaf)
    df_maf = pd.DataFrame.from_dict(dico)
    df_maf.sort_values('id', axis=0, inplace=True)
    df_maf.reset_index(drop=True, inplace=True)
    df_maf.set_index('id', drop=False, inplace=True)
    convert = np.vectorize(lambda x: convert_maf(x))
    df_maf['maf_inter'] = np.digitize(convert(df_maf['maf']),
                                      bins=np.array([0.00, 0.01, 0.05]))
    # os.remove('TMP.maf.csv')
    return df_maf


def compute_imp_err(set_names, objs, raw, idx1, idx2, sorting):
    set_err = {}
    for k, dset in dict(zip(set_names, objs)).items():
        db = para_array(raw, dset)
        print('Errors in {} dataset'.format(k).ljust(80, '.'))
        single_errors = per_element_error(db, 0)
        df = pd.DataFrame(data=single_errors, index=idx1, columns=idx2)
        set_err[k] = df
        df.to_csv(k + '{}.chunk{}.csv'.format('.sorted' if sorting else '', prm.CHK_SZ),
                  sep='\t',
                  encoding='utf-8')
    return set_err


def tag_heteroz(arr):
    if np.isin(arr, 0).any() and np.isin(arr, 1).any():
        return 1
    else:
        return 0


def tag_missing(arr):
    if np.all(arr) == -1:
        return 2
    elif np.isin(arr, -1).any():
        return 1
    else:
        return 0


def minor_allele(arr):
    if np.sum(arr) == 2:
        return 2
    elif np.isin(arr, 1).any():
        return 1
    else:
        return 0


def per_site_heteroz(vcf_path=None, gt_array=None):
    """

    :param vcf_obj:
    :return:
    """
    dic_het = collections.OrderedDict()

    if vcf_path is not None and gt_array is not None:
        raise ValueError

    elif vcf_path is None and gt_array is None:
        het = None

    elif vcf_path is not None:
        vcf_obj = VCF(vcf_path)
        for var in vcf_obj:
            gt = np.array(var.genotypes)[:, :-1]
            dic_het[var.ID] = np.sum(
                np.apply_along_axis(
                tag_heteroz, 1, gt
                ).astype(int)
            )

        het = np.asarray([[k, v] for k, v in dic_het.items()])

    else:
        tag = np.apply_along_axis(
            tag_heteroz, -1, gt_array[:, :, :-1]
        ).astype(int)
        het = np.apply_along_axis(np.sum, 0, tag)

    return het


def count_missing_alleles(vcf_path=None, gt_array=None, id='id'):
    """

    :param vcf_obj:
    :return:
    """
    dic_mis = collections.OrderedDict()

    if vcf_path is not None and gt_array is not None:
        raise ValueError

    elif vcf_path is None and gt_array is None:
        mis = None

    elif vcf_path is not None:
        vcf_obj = VCF(vcf_path)
        for var in vcf_obj:
            gt = np.array(var.genotypes)[:, :-1]
            if id == 'id':
                dic_mis[var.ID] = np.sum(
                    np.apply_along_axis(
                        tag_missing, -1, gt)
                )/(len(var.genotypes)*2)
            if id == 'chrom:pos':
                dic_mis['{}:{}'.format(var.CHROM, var.POS)] = np.sum(
                    np.apply_along_axis(
                        tag_missing, -1, gt)
                )/(len(var.genotypes)*2)
        mis = np.asarray([[k, v] for k, v in dic_mis.items()])

    else:
        tag = np.apply_along_axis(
            tag_missing, -1, gt_array[:, :, :-1]
        ).astype(int)
        mis = np.apply_along_axis(np.sum, 0, tag)

    return mis


def compute_maf(vcf_path, idt='id', verbose=False):
    """

    :param vcf_obj:
    :return:
    """
    print('Computing MAFs from {}'.format(vcf_path).ljust(80, '.'))
    vcf_obj = VCF(vcf_path)
    dic_maf = {}
    for i, var in enumerate(vcf_obj):
        if idt == 'id':
            varid = var.ID
        if idt == 'chrom:pos':
            varid = ':'.join([str(var.CHROM), str(var.POS)])
        # GT:DS:GP
        if verbose:
            try:
                print(str(i) + '\t' + var.ID + '\t' + var.format('GP'))
            except:
                print(str(i) + '\t', var.ID)
        try:
            gt = np.array(var.genotypes)[:, :-1]
            # minors = np.sum(np.apply_along_axis(minor_allele, 1, gt))
            # notmis = np.sum(~np.where(np.apply_along_axis(tag_missing, 1, gt) == 2))
            # dic_maf[varid] = minors / notmis
            # dic_maf[varid] = np.sum(np.apply_along_axis(minor_allele, 1, gt)) / (len(var.genotypes) * 2) # my version
            dic_maf[varid] = var.aaf
        except TypeError:
            pass

    return dic_maf


def compute_maf_evol(set, gtgl='gt', chk_sz=None):
    """

    :param set: str, short name of the data set pooled/missing...
    :param df_maf: maf from the original vcf-chunked file, with markers ID as index
    :return:
    """
    if chk_sz is None:
        chk_sz = prm.CHK_SZ
    print('\r\nSet --> ', set)
    # TODO: set preimp postimp path as input params
    steps = {'preimp': 'IMP.chr20.{}.snps.gt.chunk{}.vcf.gz'.format(set, str(chk_sz)),
             'postimp': 'IMP.chr20.{}.imputed.{}.chunk{}.vcf.gz'.format(set, gtgl, str(chk_sz))}
    temp = []
    for d, vcf in steps.items():
        df = pd.DataFrame.from_dict(compute_maf(vcf),
                                    orient='index',
                                    columns=[d + '_' + set])
        temp.append(df)

    return temp


def rmse_df(df, kind='rmse', ax=None, lev=None):
    df.astype(float, copy=False)
    sqr = df.applymap(lambda x: np.power(x, 2))
    if ax is None: # rmse overall
        mse = np.mean(sqr.values)
        rmse = math.sqrt(mse)
    else: # rmse per marker/per sample
        mse = sqr.mean(axis=ax, level=lev)
        rmse = mse.apply(math.sqrt).rename('rmse')
    if kind == 'rmse':
        return rmse
    else:
        return mse


def map_gt_gl(arr_in, unknown=[1/3, 1/3, 1/3]):
    gt = arr_in[:-1]
    if 0 in gt and 1 in gt:
        gl = np.array([0.0, 1.0, 0.0])
    elif 1 in gt and -1 not in gt:
        gl = np.array([0.0, 0.0, 1.0])
    elif 0 in gt:
        gl = np.array([1.0, 0.0, 0.0])
    elif -1 in gt and 1 in gt:
        gl = np.array([0.0, 0.5, 0.5])
    else: # np.nan, np.nan OR -1, -1 ?
        gl = np.array(unknown)
    return gl


def repr_gl_array(arr):
    arr2str = np.vectorize(lambda x: "%.5f" % x)
    arr_str = np.apply_along_axis(lambda x: ','.join(arr2str(x)),
                                  axis=-1,
                                  arr=arr)
    strg = np.apply_along_axis(lambda x: '\t'.join(x),
                               axis=0,
                               arr=np.asarray(arr_str.squeeze()))
    return strg


def bin_gl_converter(v_in, log=True):
    # v_in: cyvcf2 variant
    g_in = v_in.genotypes
    g_out = np.apply_along_axis(map_gt_gl, 1, g_in)
    logzero = np.vectorize(lambda x: -5.0 if x <= pow(10, -5) else math.log10(x))
    if log:
        g_out = logzero(g_out)
        # g_out = np.apply_along_axis(logzero, axis=0, arr=g_out)
    return g_out


def fmt_gl_variant(v_in, glfunc=bin_gl_converter):
    info = ';'.join([kv for kv in ['='.join([str(k), str(v)]) for k, v in v_in.INFO]])
    gl = repr_gl_array(np.array(list(map(glfunc, [v_in]))))
    toshow = np.asarray([v_in.CHROM,
                         v_in.POS,
                         v_in.ID,
                         ''.join(v_in.REF),
                         ''.join(v_in.ALT),
                         v_in.QUAL if not None else '.',
                         'PASS' if v_in.FILTER is None else v_in.FILTER,
                         info,
                         'GL',
                         gl],
                        dtype=str)
    towrite = '\t'.join(toshow) + '\n'
    return towrite


def file_likelihood_converter(f_in, f_out, func=bin_gl_converter):
    str_header = '##FORMAT=<ID=GL,Number=G,Type=Float,Description="three log10-scaled likelihoods for RR,RA,AA genotypes">'
    dic_header = {'ID': 'GL',
                  'Number': 'G',
                  'Type': 'Float',
                  'Description': 'three log10-scaled likelihoods for RR,RA,AA genotypes'}
    vcf_in = VCF(f_in)
    vcf_in.add_format_to_header(dic_header)

    for lin in vcf_in.header_iter():
        try:
            lin['ID'] == 'GT'
            del lin
        except KeyError:
            pass

    w1 = Writer(f_out, vcf_in)
    w1.write_header()
    w1.close()

    with open(f_out, 'ab') as w2:
        for var_in in vcf_in:
            stream = fmt_gl_variant(var_in, glfunc=func).encode()
            w2.write(stream)


def plot_maf_gl():
    """
    Shows how MAF = f impacts GL values for GT RR|RA|AA.
    :return:
    """
    import matplotlib.pyplot as plt

    f = np.arange(0.0, 1.0, 0.01)
    rr = np.vectorize(lambda x: (1 - x)**2)
    ra = np.vectorize(lambda x: 2 * (1 - x) * x)
    aa = np.vectorize(lambda x: x ** 2)

    plt.plot(f, rr(f), 'b-', label='RR')
    plt.plot(f, ra(f), 'g-', label='RA')
    plt.plot(f, aa(f), 'r-', label='AA')
    plt.legend()
    plt.show()


if __name__ == '__main__':
    plot_maf_gl()
