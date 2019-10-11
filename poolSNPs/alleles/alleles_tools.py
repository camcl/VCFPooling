from cyvcf2 import VCF, Writer
import numpy as np
import pandas as pd
import subprocess
import collections
from scipy.stats import *
import math

import scripts.poolSNPs.parameters as prm
from persotools.files import *

"""
Compute the imputation errors.
Vertically: read chunks of 1000 markers.
Plot different analysis.
"""


### GENERAL TOOLS

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


def vcf2dframe(vcf_obj, size, idt='id'):
    """
    Throws the genotypes values of a VCF file into an array.
    :param vcf_obj: VCF file read with cyvcf2
    :return: two dataframes, one for each allele of the genotype
    """
    vars = []
    arr = np.zeros((1, size, 2), dtype=int)
    for v in vcf_obj:
        if idt == 'id':
            vars.append(v.ID)
        if idt == 'chrom:pos':
            vars.append(':'.join([str(v.CHROM), str(v.POS)]))
        arr = np.vstack((arr,
                         np.expand_dims(np.array(v.genotypes)[:, :-1],
                                        axis=0)))
    df0 = pd.DataFrame(arr[1:,:,0], index=vars, columns=vcf_obj.samples)
    df1 = pd.DataFrame(arr[1:, :, 1], index=vars, columns=vcf_obj.samples)
    return df0, df1


def convert_aaf(x: object):
    """
    Cast AAF as float or NaN
    :param x:
    :return:
    """
    try:
        x = float(x)
    except TypeError:
        x = np.nan
    finally:
        return x


def get_aaf(vcf_raw: str, idt: str = 'id') -> pd.DataFrame:
    """
    Read the allele frequencies info field for all variants in a file.
    Return the computed alternate allele frequencies computed per variant from the input file.
    Variants must be GT filled in.
    :param vcf_raw: path to the VCF-formatted file with GT INFO fields
    :param id: kind of ID wished for pandas.Index: 'id'|'chrom:pos̈́'
    :return: id-indexed dataframe with AF from INFO-field, computed AAF and binned AF_INFO
    """
    dico = {'id': list(),
            'af_info': list(),  # alternate allele frequency estimated on the population
            'aaf': list()}  # alternate allele frequency estimated on the file
    try:
        if idt == 'id':
            for var in VCF(vcf_raw):
                dico['id'].append(var.ID)
                try:
                    dico['af_info'].append(var.INFO['AF'])
                except KeyError:
                    dico['af_info'].append(np.nan)
                dico['aaf'].append(var.aaf)
        if idt == 'chrom:pos':
            for var in VCF(vcf_raw):
                dico['id'].append(':'.join([str(var.CHROM), str(var.POS)]))
                try:
                    dico['af_info'].append(var.INFO['AF'])
                except KeyError:
                    dico['af_info'].append(np.nan)
                dico['aaf'].append(var.aaf)
        df_aaf = pd.DataFrame.from_dict(dico)
        df_aaf.sort_values('id', axis=0, inplace=True)
        df_aaf.reset_index(drop=True, inplace=True)
        df_aaf.set_index('id', drop=False, inplace=True)
        convert = np.vectorize(lambda x: convert_aaf(x))
        df_aaf['aaf_bin'] = np.digitize(convert(df_aaf['af_info']),
                                        bins=np.array([0.00, 0.01, 0.05]))
    except IOError or UnicodeDecodeError or MemoryError:
        if idt == 'id':
            subprocess.run(['''bcftools query -f '%ID\t%AF\n' {0} > {1}/TMP.aaf.csv'''.format(vcf_raw, os.getcwd())],
                           shell=True,
                           cwd=prm.DATA_PATH)
        if idt == 'chrom:pos':
            subprocess.run(['''bcftools query -f '%CHROM:%POS\t%AF\n' {0} > {1}/TMP.aaf.csv'''.format(vcf_raw, os.getcwd())],
                           shell=True,
                           cwd=prm.DATA_PATH)
        df_aaf = pd.read_csv('TMP.aaf.csv',
                             sep='\t',
                             names=['id', 'af_info'])
        os.remove('TMP.aaf.csv')

    return df_aaf


# TODO: Something wrong in that function!!!
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


def compute_discordance(setnames: list, objs: list, raw: str, idx1, idx2, sorting) -> dict:
    """
    Compute discordance per variant per sample on input data sets.
    Trinary encoding is used for GT genotypes for avoiding 4d array
    :param setnames: names for dataframes: 'pooled'/'missing'
    :param objs: dataframes to evaluate
    :param raw: dataframe with ground truth genotype values
    :param idx1: lines index
    :param idx2: columns index
    :param sorting:
    :return:
    """
    set_err = {}
    for k, dset in dict(zip(setnames, objs)).items():
        # sum the 2 alleles
        raw_trinary = raw.sum(axis=-1)
        # print('\nraw_trinary: ', describe(raw_trinary))
        dset_trinary = dset.sum(axis=-1)
        # print('\ndset_trinary: ', describe(dset_trinary))
        # compare trinary genotypes: ground truth vs. computed
        mismatch = np.subtract(raw_trinary, dset_trinary)
        print('Discordance in {} dataset'.format(k).ljust(80, '.'))
        single_errors = np.absolute(mismatch) / 2
        # print('\nsingle_errors: ', describe(single_errors))
        # save discordance result to csv
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


def per_site_zygosity(zygosity: str, vcf_path, idt: str = 'id') -> dict:
    """

    :param zygosity: 'het', 'hom_ref', 'hom_alt'
    :param vcf_path: should be GT formatted
    :param idt:
    :return:
    """
    dic_zyg = collections.OrderedDict()

    try:
        vcf_obj = VCF(vcf_path)
        if zygosity == 'het':
            for var in vcf_obj:
                dic_zyg[(var.ID, '{}:{}'.format(var.CHROM, var.POS))] = var.num_het / var.num_called
        elif zygosity == 'hom_ref':
            for var in vcf_obj:
                dic_zyg[(var.ID, '{}:{}'.format(var.CHROM, var.POS))] = var.num_hom_ref / var.num_called
        elif zygosity == 'hom_alt':
            for var in vcf_obj:
                dic_zyg[(var.ID, '{}:{}'.format(var.CHROM, var.POS))] = var.num_hom_alt / var.num_called

        if idt == 'id':
            zyg = dict([(k[0], v) for k, v in dic_zyg.items()])
        if idt == 'chrom:pos':
            zyg = dict([(k[1], v) for k, v in dic_zyg.items()])

    except IOError:
        print('VCF-file should be GT formatted')

    return zyg


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


def compute_aaf(vcf_path: FilePath, idt: str = 'id', verbose: bool = False):
    """
    Return only AAF computed on the input file, per variant ID
    :param vcf_path:
    :param idt:
    :param verbose:
    :return:
    """
    print('Computing AAFs from {}/{}'.format(os.path.dirname(vcf_path), os.path.basename(vcf_path)).ljust(80, '.'))
    print(os.getcwd())
    import inspect
    curframe = inspect.currentframe()
    print(curframe.f_code)
    print(curframe.f_back.f_locals)
    print(curframe.f_locals)
    print(curframe.f_globals)
    print('compute_aaf ID --> ', idt)
    vcf_obj = VCF(vcf_path)
    dic_aaf = {}
    for i, var in enumerate(vcf_obj):
        if idt == 'id':
            varid = var.ID
        if idt == 'chrom:pos':
            varid = ':'.join([str(var.CHROM), str(var.POS)])
        # GT:DS:GP
        if verbose:
            try:
                print(str(i) + '\t' + var.ID + '\t' + var.format('GP'))
            except TypeError:
                print(str(i) + '\t', var.ID)
        try:
            dic_aaf[varid] = var.aaf
        except TypeError:
            pass
        print(varid)

    return dic_aaf


def compute_aaf_evol(wd: str, dset: str, idt: str = 'id') -> list:
    """
    Data set for genotypes before imputation = output from phasing step
    Data set for genotypes after imputation = output from imputation step
    :param dset: short name for the data set kind e.g. pooled/missing
    :param fpre:
    :param fpost:
    :param idt:
    :return: list of DataFrames with AAF values before/after imputation
    """
    print('\r\ndset --> {} in {}'.format(dset, wd))
    os.chdir(wd)
    if dset.lower() == 'pooled':
        steps = {'preimp': os.path.join(os.path.dirname(wd), prm.POOLED['b1']) + '.vcf.gz',
                 'postimp': prm.POOLED['gtonly'] + '.vcf.gz'}
    if dset.lower() == 'missing':
        steps = {'preimp': os.path.join(os.path.dirname(wd), prm.MISSING['b1']) + '.vcf.gz',
                 'postimp': prm.MISSING['gtonly'] + '.vcf.gz'}
    temp = []
    for d, vcf in steps.items():
        df = pd.DataFrame.from_dict(compute_aaf(vcf, idt=idt),
                                    orient='index',
                                    columns=[d + '_' + dset])
        temp.append(df)

    return temp


def compute_zygosity_evol(zygosity: str, idt: str = 'id') -> list:
    """

    :param dset: str, short name of the data dset pooled/missing
    :param fpre:
    :param fpost:
    :param zygosity:
    :param idt:
    :return: list of DataFrames with het values for each variant before/after imputation
    """
    print('\r\nZygosity --> ', zygosity)
    if zygosity.lower() == 'pooled':
        steps = {'preimp': prm.POOLED['b1'] + '.vcf.gz',
                 'postimp': prm.POOLED['gtonly'] + '.vcf.gz'}
    if zygosity.lower() == 'missing':
        steps = {'preimp': prm.MISSING['b1'] + '.vcf.gz',
                 'postimp': prm.MISSING['gtonly'] + '.vcf.gz'}
    temp = []
    for d, vcf in steps.items():
        df = pd.DataFrame.from_dict(per_site_zygosity(zygosity, vcf, idt=idt),
                                    orient='index',
                                    columns=[d + '_' + zygosity])
        df.set_index(idt, inplace=True)
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


def extract_variant_onfreq(file_in, aaf_range):
    """
    Returns variants where the alternate allele has a frequency
    in the input range.
    :param file_in: VCF-formatted file
    :param aaf_range: list: [min, max] of the range the AAF should be comprised in.
    :return: dataframe of cyvcf2.Variants
    """
    aafs = get_aaf(file_in, id='chrom:pos')
    print(aafs.describe())
    inf, sup = aaf_range[0], aaf_range[1]
    #extracted = aafs.loc[(aafs['af_info'] >= min and aafs['af_info'] <= max).all()]
    extracted = aafs.query('af_info >= @inf & af_info <= @sup')

    return extracted


def map_gt_gl(arr_in, unknown=[1/3, 1/3, 1/3]):
    gt = arr_in[:-1]
    if 0 in gt and 1 in gt:
        gl = np.array([0.0, 1.0, 0.0])
    elif 1 in gt and -1 not in gt:
        gl = np.array([0.0, 0.0, 1.0])
    elif 0 in gt and -1 not in gt:
        gl = np.array([1.0, 0.0, 0.0])
    elif -1 in gt and 1 in gt:
        gl = np.array([0.0, 0.5, 0.5])
    elif -1 in gt and 0 in gt:
        gl = np.array([0.5, 0.5, 0.0])
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


def get_pop():
    df_pop = pd.read_csv(os.path.join(prm.DATA_PATH, '20130606_sample_info.csv'),
                         sep=',',
                         header=0,
                         usecols=[0, 1])
    df_pop.sort_values(['Sample', 'Population'], axis=0, inplace=True)
    return df_pop


def make_index(raw_data, src=None, idt='id'):
    if src is None:
        src = os.path.join(prm.WD, 'gt', 'stratified', prm.CHKFILE)
    group = get_pop()
    df_aaf = get_aaf(src, id=idt).loc[:, ['id', 'af_info', 'aaf_bin']]

    samples = pd.Series(VCF(raw_data).samples,
                        name='Sample').str.rstrip('_IMP').to_frame()
    df_pop = group.merge(samples, on='Sample', how='inner')

    aaf_idx = pd.MultiIndex.from_arrays(list(zip(*df_aaf.values)), names=['id', 'af_info', 'aaf_bin'])
    pop_idx = pd.MultiIndex.from_arrays(list(zip(*df_pop.values)), names=['Sample', 'Population'])

    return aaf_idx, pop_idx


def normalize(v):
    norm = np.linalg.norm(v)
    if norm == 0:
        return v
    return v / norm



