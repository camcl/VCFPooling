import subprocess
import warnings
import numpy as np
import pandas as pd
import time
import math
from scipy.stats import bernoulli as bn
from cyvcf2 import VCF, Writer, Variant

from scripts.VCFPooling.poolSNPs import parameters as prm
from scripts.VCFPooling.poolSNPs.alleles import alleles_tools as alltls
from scripts.VCFPooling.python.archived import frequency_distribution as allfqc
from scripts.VCFPooling.poolSNPs import pybcf

from persotools.files import *
from typing import *

global SAMPLES

warnings.simplefilter('ignore')

nb_cores = os.cpu_count()

# VCF file to be pooled
data = VCF(os.path.join(prm.WD, 'gt', prm.CHKFILE), threads=nb_cores)

SAMPLES = data.samples

### README

"""
README: 
Requirements: 
- cyvcf2
- bcftools

Usage for creating pooled genotypes likelihoods
from a set of samples in a VCF file:
(1) Instanciate a Pool structure:
    pool = Pool()
(2) Update the pool's parameters with file's data:
    pool.set_subset(...)
(3) Fill the pooling matrix with data from file:
    pool.set_line_values(...)
(4) Implement the needed class methods
    ex. pool.decode_genotypes_gt(...)
"""

# ## TOOLS FUNCTIONS


def random_delete(rate: float = 0.01, activate: bool = False) -> bool:
    """
    Delete randomly genotypes in a variant call.
    :param arr: variant.genotypes instance
    :param rate: proportion of genoyptes to delete over the call
    :param activate: activation of the deletion function
    :return: array of boolean indicating if the genotypes have to be deleted (1) or not (0).
    """
    if activate:
        flag: bool = bn.rvs(p=rate, size=1)
    else:
        flag: bool = False

    return flag


# ## CLASS AND METHODS FOR POOLING A SET OF SAMPLES


def split_pools(idv_id: List[str], pool_size: int, seed: object = None) -> Tuple[List[str], List[str]]:
    """
    Assign samples to pools and output the remaning samples that could not be assigned (to few samples).
    :param idv_id: sample's identiifer in the VCF-file.
    :param pool_size: size of group expected (e.g. 16 if 4*4 pooling design)
    :return: samples identifiers grouped and remaining samples
    """
    idv_nb = len(idv_id)
    nb_pool, left = divmod(idv_nb, pool_size)
    print('Number of {}-sized pools to create, number of single samples remaining: '.format(pool_size),
          nb_pool, left)
    # Create random pools
    if seed is not None:
        np.random.seed(seed)
    pools = np.random.choice(idv_id, size=(nb_pool, pool_size), replace=False)
    left = np.isin(idv_id, pools.flatten())
    singles = np.extract(~left, idv_id)
    pools = pools.tolist()

    return pools, singles


class SNPsPool(np.ndarray):
    """
    Simulates the different steps of a genotype pooling process.
    Builds the pooling design matrix.
    Encodes pooled genotypes per group.
    Decodes encoded genotypes of the pools back to individual genotypes.
    """
    def __new__(cls,
                shape: tuple = (4, 4),
                id_len: int = 8,
                pools_nb: int = 8,
                pools_size: int = 4) -> np.ndarray:
        """
        Define the basic structure for a pool i.e.
        a squared matrix to fill with the variables IDs/GT/GL.
        :param shape: tuple, shape of the pool
        :param id_len: max number of char of the variables IDs
        :return: a matrix with dims 'shape', filled with str types
        """
        cls.id_len = id_len
        id = 'U' + str(cls.id_len)
        cls.pools_nb = pools_nb
        cls.pools_size = pools_size
        return np.empty_like(super(SNPsPool, cls).__new__(cls, shape),
                             dtype=id)

    def design_matrix(self, random: bool = False) -> np.ndarray:
        """
        That function is not intended to be called explicitly.
        :param random: bool for dispatching idv randomly in the matrix?
        :return: design matrix. Numpy array.
        """
        pools_size: int = self.pools_size
        design: np.ndarray = np.zeros((self.pools_nb, self.size), dtype=int)
        if not random:
            for i in range(int(self.pools_nb/self.ndim)):
                j = i * pools_size
                design[i, j:j+pools_size] = [1]*pools_size
            for i in range(int(self.pools_nb/self.ndim), self.pools_nb):
                j = i - pools_size
                design[i,
                       [j+k*pools_size for k in range(pools_size)]] = 1
        return design

    def set_subset(self, subset: np.ndarray) -> np.ndarray:
        # TODO: refactor with more sel-explanatory name
        """
        Fills the pooling matrix according to the (default) design
        and the input list of samples.
        :param subset: 1D-nparray-like object with the variables IDs
        :return: pooling matrix with samples' names.
        """
        self.__setattr__('subset', subset)  # from split_pools
        sub = self.__getattribute__('subset')
        try:
            for i in range(self.shape[0]):
                self[i, :] = sub[:self.shape[1]]
                sub = sub[self.shape[1]:]
        except Exception as exc:
            if len(self.subset) > self.size:
                raise ValueError('The input you gave is too long') from exc
            if len(self.subset) < self.size:
                raise ValueError('The input you gave is too short') from exc
            if type(self.subset) != np.ndarray and type(self.subset) != list:
                raise TypeError('The input is not a 1D-array-like') from exc
            if len(self.subset) > 0 and type(self.subset[0]) != str:
                raise TypeError('The input does not contain str-type elements') from exc

        return self

    def get_subset(self) -> np.ndarray:
        """
        Flatten the matrix of pooled samples identifiers.
        :return: flattened array of samples identifiers or genotypes.
        """
        ids = self.flatten()  # .reshape(1, self.size)
        return ids

    def pools_list(self) -> List[str]:
        """
        Put samples into matrix structure, just for representation
        :return: list of matrix-shaped samples identifiers
        """
        design = self.design_matrix()
        if np.where(self == '', False, True).all():
            pools_: List[str] = []
            for i in range(design.shape[0]):
                cache = (design[i, :].reshape(self.shape) == False)  # cache = ~design[i, :].reshape(self.shape)
                pool = np.ma.masked_array(self, mask=cache)
                pools_.append(pool.compressed())
            return pools_
            # return just for being able to print the list if wished

    def set_line_values(self, samples: list, variant: Variant,
                        sig: object = None,
                        params: List[float] = [], interp: object = None) -> None:
        # TODO: refactor with property decorator
        """
        Attach sigmoid-transformed alternate allele frequencies to the current variant.
        :param samples: samples identifiers from the VCF-file
        :param variant: variant identifier from the VCF-file
        :param sig: sigmoid object computed from another set of pooled genotype data
        :param params: parameters to pass for the sigmoid approximation
        :param interp: interpolated sigmoid object for sigmoid approximation
        :return: variant object with attributes attached, e.g. twisted alternate allele frequency
        """
        self.__setattr__('variant', variant.genotypes)
        self.__setattr__('samples', samples)
        self.__setattr__('var_pos', str(variant.POS))
        self.__setattr__('aaf', variant.aaf)
        if sig is not None:
            self.__setattr__('aat', sig.call_sigmoid(params, self.aaf))
            self.__setattr__('aat_', sig.call_sigmoid_derivative(interp, self.aaf))
        else:
            self.__setattr__('aat', np.nan)
            self.__setattr__('aat_', np.nan)

    def get_call(self) -> np.ndarray:
        """
        Get pooled samples genotypes for a given variant
        :return:
        """
        subs = self.get_subset()
        idx = np.argwhere(np.isin(self.samples, subs))
        self.__setattr__('call', np.asarray(self.variant)[idx])
        return self.call

    def pool_genotypes(self) -> np.ndarray:
        """
        Computes genotypes of the different pools.
        :return: array of {0,1} GLs for RR|RA|AA for each pool
        """
        call: np.ndarray = self.get_call().reshape((1, self.size, 3))
        scores: np.ndarray = np.apply_along_axis(sum, axis=-1, arr=call[:, :, :-1])

        if np.isin(call, -1).any():
            x = np.ones((1, self.pools_nb, 1))
            y = np.asarray([-1, -1, 0])
            b = np.broadcast(x, y)
            p = np.empty(b.shape)
            p.flat = [u * v for (u, v) in b]
        else:
            pooled_gt = np.dot(self.design_matrix(),
                               np.transpose(scores)).reshape((1, self.pools_nb, 1))
            pooled_gt = np.broadcast_to(pooled_gt, (1, self.pools_nb, 3))
            p = np.apply_along_axis(self.pooler_gt, axis=-1, arr=pooled_gt)

        return p  # list of gt for the 8 pools from design matrix

    def pooler_gt(self, a: np.ndarray) -> np.ndarray:
        """
        Decodes pooled scores into individual GT.
        :param a: score from matrix-vector pooling
        :return: pool's true genotype with phase
        """
        if np.all(a == 0):  # RR * RR * RR * RR
            gt = [0, 0, 0]
        elif np.all(a == self.pools_nb):  # AA * AA * AA * AA
            gt = [1, 1, 0]
        else:
            gt = [1, 0, 0]
        return gt

    def decode_genotypes_gt(self, samples_gt: np.ndarray) -> np.ndarray:
        """
        Recomputes true genotypes of samples with/without pooling/missing data
        :param pooled_samples: Variant.genotypes (unpooled samples' true genotypes with phase)
        :return: individual samples genotypes (true genotype with phase)
        """
        pooled: np.ndarray = self.pool_genotypes()  # pooled[:, :, -1]: bool = phase of the genotype
        scores: np.ndarray = np.apply_along_axis(sum, axis=-1, arr=pooled[:, :, :-1])
        p = np.argwhere(np.isin(self.samples, self.subset))

        count_alt: Callable[int, int] = lambda x: 1 if 1 in x else 0
        count_ref: Callable[int, int] = lambda x: 1 if 0 in x else 0

        alt_row: int = np.sum(np.apply_along_axis(count_alt,
                                                  axis=-1,
                                                  arr=pooled[:, :4, :-1]))
        alt_col: int = np.sum(np.apply_along_axis(count_alt,
                                                  axis=-1,
                                                  arr=pooled[:, 4:, :-1]))
        ref_row: int = np.sum(np.apply_along_axis(count_ref,
                                                  axis=-1,
                                                  arr=pooled[:, :4, :-1]))
        ref_col: int = np.sum(np.apply_along_axis(count_ref,
                                                  axis=-1,
                                                  arr=pooled[:, 4:, :-1]))

        nb_alt: int = alt_row + alt_col
        nb_ref: int = ref_row + ref_col

        if np.isin(pooled, -1).any():
            x = np.ones((1, self.size, 1))
            y = np.asarray([-1, -1, 0])
            b = np.broadcast(x, y)
            decoded_gt = np.empty(b.shape)
            decoded_gt.flat = [u * v for (u, v) in b]

        else:
            encoded = np.dot(scores,
                             self.design_matrix()).reshape(1, self.size, 1)
            b = np.broadcast_to(encoded, (1, self.size, 3))
            if nb_alt == 0:
                decoded_gt = np.zeros_like(b)
            elif nb_ref == 0:
                aa = np.array([1, 1, 0])
                decoded_gt = np.tile(aa, self.size).reshape((1, self.size, 3))
            elif nb_alt == 2:
                decoder: Callable = lambda x: [1, -1, 0] if np.all(x == 2) else [0, 0, 0]
                # np.all() because of b.shape
                decoded_gt = np.apply_along_axis(decoder, axis=-1, arr=b)
            elif nb_ref == 2:  # symmetric case with ref allele in only 2 pools: individual is RR or RA
                decoder: Callable = lambda x: [0, -1, 0] if np.all(x == 2) else [1, 1, 0]
                decoded_gt = np.apply_along_axis(decoder, axis=-1, arr=b)
            else:  # nb_alt > 2 and nb_alt < 8: # nb_alt = 2*n with n > 1
                decoded_gt = np.apply_along_axis(self.multidecoder_gt, axis=-1, arr=b)

        np.put_along_axis(samples_gt,
                          np.broadcast_to(p, (self.size, 3)),
                          decoded_gt.squeeze(),
                          axis=0)

        return samples_gt

    def decode_genotypes_gl(self, samples_gl: np.ndarray, dict_gl: dict) -> np.ndarray:
        """
        Recomputes genotypes likelihoods of samples with/without pooling/missing data
        :param samples_gl: samples' true genotypes with phase
        :param dict_gl: likelihoods values to set when encountering missing genotypes
        :return: individual samples genotypes (genotype likelihoods)
        """
        samples_gl = samples_gl.astype(float)  # avoid truncating GL
        pooled: np.ndarray = self.pool_genotypes()  # outputs unphased genotypes
        scores: np.ndarray = np.apply_along_axis(sum, axis=-1, arr=pooled[:, :, :-1]).flatten()
        p = np.argwhere(np.isin(self.samples, self.subset))
        count_alt: Callable[int] = lambda x: 1 if 1 in x else 0
        count_ref: Callable[int] = lambda x: 1 if 0 in x else 0

        alt_row: int = np.sum(np.apply_along_axis(count_alt,
                                                  axis=-1,
                                                  arr=pooled[:, :4, :-1]))  # :-1 excludes the phase boolean
        alt_col: int = np.sum(np.apply_along_axis(count_alt,
                                                  axis=-1,
                                                  arr=pooled[:, 4:, :-1]))
        ref_row: int = np.sum(np.apply_along_axis(count_ref,
                                                  axis=-1,
                                                  arr=pooled[:, :4, :-1]))
        ref_col: int = np.sum(np.apply_along_axis(count_ref,
                                                  axis=-1,
                                                  arr=pooled[:, 4:, :-1]))

        rowcounts, colcounts = self.rowcolcounts(scores)

        nb_alt: int = alt_row + alt_col
        nb_ref: int = ref_row + ref_col

        if np.isin(pooled, -1).any():
            x = np.ones((1, self.size, 1))
            y = np.asarray([1/3, 1/3, 1/3])
            b = np.broadcast(x, y)
            decoded_gl = np.empty(b.shape)
            decoded_gl.flat = [u * v for (u, v) in b]

        else:
            encoded = np.dot(scores,
                             self.design_matrix()).reshape(1, self.size, 1)
            b = np.broadcast_to(encoded, (1, self.size, 3)).astype(float)
            if nb_alt == 0:
                rr = np.array([1, 0, 0])
                decoded_gl = np.tile(rr, self.size).reshape((1, self.size, 3))
            elif nb_ref == 0:
                aa = np.array([0, 0, 1])
                decoded_gl = np.tile(aa, self.size).reshape((1, self.size, 3))
            else:  # nb_alt >= 2 and nb_alt < 8: # nb_alt = 2*n with n > 1
                unknown = dict_gl[tuple([*rowcounts, *colcounts, 1, 1])]  # 1,1 only cross-section that outputs unknown genotypes
                decoded_gl = np.apply_along_axis(self.multidecoder_gl, -1, b, unknown)

        np.put_along_axis(samples_gl,
                          np.broadcast_to(p, (self.size, 3)),
                          decoded_gl.squeeze(),
                          axis=0)

        return samples_gl

    @staticmethod
    def rowcolcounts(a: np.ndarray) -> Tuple[Tuple[int, int, int]]:
        """
        Count number of pooled RR|RA|AA genotypes over all rows and columns of a pooled matrix
        :param a: score i.e. trinary-encoded true genotypes
        :return: counts of genotypes for the rows and columns
        """
        # a should be scores
        count_rr: Callable[int] = np.vectorize(lambda x: 1 if x == 0 else 0)
        count_aa: Callable[int] = np.vectorize(lambda x: 1 if x == 2 else 0)
        count_ra: Callable[int] = np.vectorize(lambda x: 1 if x == 1 else 0)

        rowcounts: List[int] = [0, 0, 0]  # counts number of RR, RA, AA rows
        colcounts: List[int] = [0, 0, 0]  # counts number of RR, RA, AA columns
        rowcounts[0] = np.sum(count_rr(a[:4]))
        colcounts[0] = np.sum(count_rr(a[4:]))
        rowcounts[1] = np.sum(count_ra(a[:4]))
        colcounts[1] = np.sum(count_ra(a[4:]))
        rowcounts[2] = np.sum(count_aa(a[:4]))
        colcounts[2] = np.sum(count_aa(a[4:]))

        return tuple(rowcounts), tuple(colcounts)

    @staticmethod
    def multidecoder_gt(a: np.ndarray) -> prm.GTtype:
        """
        Decodes pooled scores into individual GT.
        :param a: score
        :return: true genotype with phase
        """
        if np.all(a == 2):  # RA * RA
            gt = [-1, -1, 0]
        elif np.all(a == 1) or np.all(a == 0):  # RA * RR or RR * RR
            gt = [0, 0, 0]
        else:
            gt = [1, 1, 0]
        return gt

    @staticmethod
    def multidecoder_gl(a: np.ndarray, mis: List[float]) -> prm.GLtype:
        """
        Decodes pooled scores into individual GL.
        :param a: score
        :param mis: GL to set when unknown genotype
        :return: genotypes likeihoods
        """
        if np.all(a == 2):  # RA * RA
            gl = mis
        elif np.all(a == 1) or np.all(a == 0):  # RA * RR or RR * RR
            gl = [1.0, 0.0, 0.0]
        else:
            gl = [0.0, 0.0, 1.0]
        return gl

    @staticmethod
    def normalize(v):
        """
        Normalize (Z-norm) an array. See NumPy documentation.
        :param v: array to normalize
        :return: normalized array
        """
        norm = np.linalg.norm(v)
        if norm == 0:
            return v
        return v / norm

    def twister(self, a: prm.GLtype, f0: float, f1: float, f1_: float) -> prm.GLtype:
        """
        Adjusts the GLs according to the AAF gap due to pooling
        :param a: genotype likelihoods
        :param f0: pooled allele frequency value
        :param f1: sigmoid reciprocal of the pooled allele frequency
        :param f1_: derivative value
        :return: genotypes likelihoods
        """
        if f0 >= 0.4 and f0 < 0.6:
            t = np.multiply(a, [1, max((f1/f0), (f0/f1))**f1_, 1])
        elif f0 >= 0.6 and f1 > 0.001:
            t = np.multiply(a, [math.fabs(f1-f0)*10, (f1/f0)**f1_, (f0/f1)**f1_])
        elif f0 > 0.001 and f0 < 0.4 and f1 > 0.001:
            t = np.multiply(a, [(f0/f1), (f0/f1)**f1_, math.fabs(f1-f0)*10])
        else:
            t = np.multiply(a, [(1 + f0), 1, f0])

        return self.normalize(t)

    def __array_finalize__(self, obj: object) -> None:
        """
        Constructor needed for subclassing NumPy arrays.
        See online documentation.
        :param obj:
        :return:
        """
        if obj is None: return
        self.info = getattr(obj, 'info', None)


# ## PROCESS SUMMARY


def process_file(data: VCF, groups: list, simul: str) -> None:
    #TODO: clean/refactor execution comments like processed file name
    """
    Computes and rewrites genotypes of all individuals for all samples from input files
    :param data: cyvcf2 object reader pointing on a VCF-file
    :param groups: samples identifiers split in pools
    :param f: integer, index of the file to process in the list
    :param fileout: VCF-files with simulated pooled or randomly missing genotypes
    """
    print('Simulation type: ', simul)
    print('File out: ', os.path.join(os.getcwd(), prm.PATH_OUT[simul]))  # prm.PATH_OUT[simul]
    if prm.GTGL == 'GL' and prm.unknown_gl == 'adaptive':
        dic_header = {'ID': 'GL',
                      'Number': 'G',
                      'Type': 'Float',
                      'Description': 'three log10-scaled likelihoods for RR,RA,AA genotypes'}
        data.add_format_to_header(dic_header)
        whead = Writer(prm.PATH_OUT[simul], data)
        whead.write_header()
        whead.close()
        w = open(prm.PATH_OUT[simul], 'ab')
        # Load adaptive GL values for missing data
        df = pd.read_csv(os.path.join(prm.WD, 'adaptive_gls.csv'),
                         header=None,
                         names=['rowsrr', 'rowsra', 'rowsaa', 'colsrr', 'colsra', 'colsaa',
                                'n', 'm',
                                'rr', 'ra', 'aa']
                         )
        df2dict = dict(((int(rwrr), int(rwra), int(rwaa), int(clrr), int(clra), int(claa),
                         int(n), int(m)),
                        [rr, ra, aa]) for rwrr, rwra, rwaa, clrr, clra, claa,
                                          n, m,
                                          rr, ra, aa in df.itertuples(index=False, name=None))

        sig = allfqc.SigmoidInterpolator(os.path.join(prm.PATH_GT_FILES, prm.RAW['gz'].replace('gl', 'gt')),
                                         os.path.join(prm.PATH_GT_FILES, prm.POOLED['gz'].replace('gl', 'gt')))
        params = sig.get_sigmoid_params()
        interp = sig.interpolate_derivative()

    else:  # prm.GTGL == 'GT' or fixed GL
        w = Writer(prm.PATH_OUT[simul], data)
        w.set_threads(4)
        df2dict = None
        sig = None
        params = None
        interp = None

    tm = time.time()
    # for n, variant in enumerate(data('20:59973567-59973568')):
    for n, variant in enumerate(data):
        process_line(groups, simul, w, variant, df2dict, sig, params, interp)
        if n % 1000 == 0:
            print('{} variants processed in {:06.2f} sec'.format(n+1, time.time()-tm).ljust(80, '.'))
        # if n+1 == 1000:
        #     break
    w.close()

    # GL converted from GT, missing GLs will be filled with [0.33, 0.33, 0.33]
    if prm.GTGL == 'GL' and prm.unknown_gl != 'adaptive':
        alltls.file_likelihood_converter(os.path.join(prm.PATH_GT_FILES,
                                                      prm.PATH_OUT[simul].replace('.gl', '.gt')) + '.gz',
                                         prm.PATH_OUT[simul])


def process_line(groups: list, simul: str, w: Writer, v: Variant, dict_gl: dict,
                 sig: object, params: List[float], interp: object, write: bool = True) -> None:
    #TODO: comments to add
    """
    From currently pointed variant object:
    Computes and rewrites genotypes of all individual for one sample.
    :param f: integer, index of the file to process in the list
    :param v: cyvcf2.cyvcf2.Variant object
    :param w: cyvcf2.cyvcf2.Writer object
    :return: variant object with pooled values for GT/GL
    """
    var = v # copy the variant object to make it readable several times
    pooled_samples = np.asarray(var.genotypes)
    sets = []
    for gp in groups[0]:
        sets.append(SNPsPool().set_subset(gp))

    if simul == 'pooled':  # sig might be not None if adaptive GL
        i = 1
        for p in sets:
            i += 1
            p.set_line_values(SAMPLES, var, sig, params, interp)
            if prm.GTGL == 'GL' and prm.unknown_gl == 'adaptive':
                    pooled_samples = p.decode_genotypes_gl(pooled_samples,
                                                           dict_gl)
            else: # prm.GTGL == 'GT' or fixed GL
                pooled_samples = p.decode_genotypes_gt(pooled_samples)

    else:  # randomly missing simulation
        for p in sets:
            p.set_line_values(SAMPLES, var)
            dlt = random_delete(activate=True)
            idx = np.argwhere(np.isin(SAMPLES, p))
            if dlt:
                if prm.GTGL == 'GL' and prm.unknown_gl == 'adaptive':
                    pooled_samples = pooled_samples.astype(float)  # avoid truncating GL
                    np.put(pooled_samples, idx, np.asarray([1/3, 1/3, 1/3]))
                else:
                    np.put(pooled_samples, idx, np.asarray([-1, -1, 0]))

    if write:
        if prm.GTGL == 'GL' and prm.unknown_gl == 'adaptive':
            # cyvcf2.Variant.genotypes does not handle GL-format
            # customize line format for GL
            logzero = np.vectorize(lambda x: -5.0 if x <= pow(10, -5) else math.log10(x))
            info = ';'.join([kv for kv in ['='.join([str(k), str(v)]) for k, v in var.INFO]])
            gl = alltls.repr_gl_array(logzero(pooled_samples))
            toshow = np.asarray([var.CHROM,
                                 var.POS,
                                 var.ID,
                                 ''.join(var.REF),
                                 ''.join(var.ALT),
                                 var.QUAL if not None else '.',
                                 'PASS' if var.FILTER is None else var.FILTER,
                                 info,
                                 'GL',
                                 gl],
                                dtype=str)
            towrite = '\t'.join(toshow) + '\n'
            stream = towrite.encode()
            w.write(stream)
        else:
            # cyvcf2.Variant.genotypes does handle GT-format
            var.genotypes = pooled_samples.tolist()
            w.write_record(var)


def init_chunk(WD: str, path_in: str, chunk: bool = True, strat: bool = False) -> Tuple[List[str], List[str]]:
    """
    Randomly select a limited amount of markers from an input VCF file.
    Obs! Long to run when stratifying markers (sorting/filtering tasks involved)
    :param chunk: boolean. If True, a new subfile from parameters.CHK_SZ
    randomly drawm markers is created i.e. a new set of SNPs.
    :return: list of 16-sized groups of samples i.e. pools
    """
    raw = VCF(os.path.join(prm.WD, 'gt', path_in), threads=nb_cores)  # VCF iterator object
    splits = split_pools(raw.samples, 16, seed=123)  # list of lists
    with open(os.path.join(prm.WD, 'gt', 'ALL.chr20.snps.allID.txt'), 'w') as f:
        for s in splits[0]:
            for i in s:
                f.write(i + os.linesep)

    if chunk:
        if not strat:
            delete_file('ALL.chr20.snps.gt.chunk{}.vcf.gz'.format(str(prm.CHK_SZ)))
            delete_file('chunk_{}.vcf'.format(str(prm.CHK_SZ)))
            delete_file('TMP.chr20.snps.gt.chunk.vcf')

            print('Extracting header'.ljust(80, '.'))
            pybcf.extract_header('ALL.chr20.snps.gt.vcf.gz', 'headers.ALL.chr20.snps.gt.vcf', WD)

            print('Sampling markers'.ljust(80, '.'))
            pybcf.chunk_markers(prm.SRCFILE, prm.CHK_SZ, WD)
            pybcf.concatenate(['headers.ALL.chr20.snps.gt.vcf', 'chunk_{}.vcf '.format(str(prm.CHK_SZ))],
                              'TMP.chr20.snps.gt.chunk.vcf', WD)

            print('BGzipping chunk file'.ljust(80, '.'))
            pybcf.bgzip('TMP.chr20.snps.gt.chunk.vcf', prm.CHKFILE, WD)
            delete_file('chunk_{}.vcf'.format(str(prm.CHK_SZ)))
            delete_file('TMP.chr20.snps.gt.chunk.vcf')
            delete_file('headers.ALL.chr20.snps.gt.vcf'.format(str(prm.CHK_SZ)))

        else:
            delete_file('ALL.chr20.snps.gt.chunk{}.vcf.gz'.format(str(prm.CHK_SZ)))
            mkdir(os.path.join(prm.DATA_PATH,
                               'gt',
                               'stratified'))
            pybcf.stratified_aaf_sampling(prm.SRCFILE,
                                          os.path.join(prm.DATA_PATH,
                                                       'gt',
                                                       'stratified'))


        # Eliminate the remaining samples (not involved in any pool)
        os.chdir(prm.PATH_GT_FILES)
        subprocess.run(' '.join(['bcftools view -Ov -o ALL.chr20.snps.gt.chunk{}.vcf'.format(str(prm.prm.CHK_SZ)),
                                 '-s',
                                 '^' + ','.join(splits[1]),
                                 'ALL.chr20.snps.gt.chunk{}.vcf.gz'.format(str(prm.CHK_SZ))]),
                       shell=True,
                       cwd=prm.PATH_GT_FILES)
        # back to vcf un-bgzipped for avoiding "file truncated" error
        pybcf.bgzip('ALL.chr20.snps.gt.chunk{}.vcf'.format(str(prm.CHK_SZ)), prm.CHKFILE, WD)
        pybcf.sort(prm.CHKFILE, WD)
        pybcf.index(prm.CHKFILE, WD)
        delete_file('ALL.chr20.snps.gt.chunk{}.vcf'.format(str(prm.CHK_SZ)))

    # Reorder samples according to pools (reorder columns int the VCF file)
    os.chdir(prm.PATH_GT_FILES)
    pybcf.bgzip(prm.CHKFILE,
                'ALL.chr20.snps.gt.chunk{}.unordered.samples.vcf.gz'.format(str(prm.CHK_SZ)),
                prm.PATH_GT_FILES)
    pybcf.index(prm.CHKFILE, prm.PATH_GT_FILES)
    pybcf.sampling('ALL.chr20.snps.gt.chunk{}.unordered.samples.vcf.gz'.format(str(prm.CHK_SZ)),
                   prm.CHKFILE,
                   '{}/ALL.chr20.snps.allID.txt'.format(prm.WD + '/gt/'),
                   prm.PATH_GT_FILES)
    pybcf.index(prm.CHKFILE, prm.PATH_GT_FILES)
    delete_file('ALL.chr20.snps.gt.chunk{}.unordered.samples.vcf.gz'.format(str(prm.CHK_SZ)))
    delete_file('ALL.chr20.snps.gt.chunk{}.unordered.samples.vcf.gz.csi'.format(str(prm.CHK_SZ)))
    os.chdir(WD)

    return splits


def run(splits: list, sim: str) -> None:
    #TODO: that function should take as input params:
    # * type of simulation for the missing data: replace iteration + modify comments,
    # * output file path,
    # *
    """

    :param splits:
    :param sim:
    :return:
    """
    vcf = VCF(os.path.join(prm.PATH_GT_FILES, prm.CHKFILE), threads=nb_cores)
    print('File to write to --> ', prm.PATH_OUT[sim])
    start = time.time()
    process_file(vcf, splits, sim)
    stop = time.time()
    print('Elapsed time for pooling the VCF files: {:06.2f} sec'.format(stop - start))


def write_truncate_vcf(path_in: str, path_out: str, trunc: int) -> int:
    """

    :param path_in:
    :param path_out:
    :param trunc:
    :return:
    """
    w = Writer(path_out, VCF(path_in, threads=nb_cores))
    for i, v in enumerate(VCF(path_in, threads=nb_cores)):
        if i == trunc:
            break
        else:
            w.write_record(v)
    return i


def subset_chunked_vcf(wd: str, src: str, path_out: str, chz_sz: int, trc: int) -> None:
    """
    
    :param wd:
    :param src:
    :param path_out:
    :param chz_sz:
    :param trc:
    :return:
    """
    for fi in path_out:
        delete_file(fi.replace('chunk' + str(chz_sz), 'chunk' + str(trc)) + '.gz.csi')
        print('Creating subchunk file for {}'.format(fi).ljust(80, '.'))
        lgth = write_truncate_vcf(fi, fi.replace('chunk' + str(prm.CHK_SZ), 'chunk' + str(trc)), trc)
        print('Subchunk size: ', lgth)

    delete_file(src.replace('chunk' + str(chz_sz), 'chunk' + str(trc)) + '.csi')
    delete_file(src[:-3].replace('chunk' + str(chz_sz), 'chunk' + str(trc)))
    print('Creating subchunk file for {}'.format(src).ljust(80, '.'))
    lgth = write_truncate_vcf(src, src[:-3].replace('chunk' + str(chz_sz), 'chunk' + str(trc)), trc)
    print('Subchunk size: ', lgth)
    # source[:-3] deleting '.gz' in the name

    print('BGzipping subchunk file from {}'.format(src).ljust(80, '.'))
    subprocess.run('bcftools view -Oz -o ALL.chr20.snps.gt.chunk{}.vcf.gz '.format(str(trc))
                   + 'ALL.chr20.snps.gt.chunk{}.vcf '.format(str(trc)),
                   shell=True,
                   cwd=wd)
    print('Sorting and indexing subchunk file'.ljust(80, '.'))
    subprocess.run('bcftools sort -Oz -o ALL.chr20.snps.gt.chunk{}.vcf.gz '.format(str(trc))
                   + 'ALL.chr20.snps.gt.chunk{}.vcf.gz'.format(str(trc)),
                   shell=True,
                   cwd=wd)
    subprocess.run('bcftools index -f ALL.chr20.snps.gt.chunk{}.vcf.gz'.format(str(trc)),
                   shell=True,
                   cwd=wd)
