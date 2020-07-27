import os
import numpy as np
cimport numpy as np
import cython
from src.VCFPooling.poolSNPs.pooler import Design
import pandas as pd


"""
Try to optimize decoding with NumPy style classes
@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
@cython.boundscheck(False)
@cython.wraparound(False)
"""

DTYPE = np.float
ctypedef double DTYPE_t
ctypedef np.npy_intp SIZE_t


class SingleBlockDecoder(object):
    """
        Simulate edeoding step in pooling.
        Proceed block-wise.
        This version is based on the use of dot-products of NumPy arrays for
         representing the lookup table with the adaptive GP values.
        """

    def __init__(self, design1block, lookup_keys, lookup_vals, fmt):
        self.ds = design1block  # block-wise
        self.fmt = fmt.upper()
        assert (self.fmt == 'GT' or self.fmt == 'GP'), 'Pooling to other formats than GT or GP not implemented'
        #cdef np.ndarray[np.int_t, ndim=2] D
        #cdef double[:, :] V
        self.D = lookup_keys
        self.V = lookup_vals

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def decode_genotypes_gt(self, pooled):
        """
        Recomputes true genotypes of samples with/without pooling/missing data
        :param pooled: sum of alleles of pooled true genotypes (unpohased)
        :return: individual samples genotypes (true genotype unphased)
        """
        count_alt = np.where(pooled >= 1, 1, 0)
        count_ref = np.where(pooled <= 1, 1, 0)

        alt_row = int(np.sum(count_alt[4:]))
        alt_col = int(np.sum(count_alt[4:]))
        ref_row = int(np.sum(count_ref[:4]))
        ref_col = int(np.sum(count_ref[4:]))

        cdef int nb_alt = alt_row + alt_col
        cdef int nb_ref = ref_row + ref_col

        encoded = np.dot(pooled, self.ds).reshape(1, self.ds.shape[1], 1)
        b = np.broadcast_to(encoded, (1, self.ds.shape[1], 2))
        if nb_alt == 0:
            decoded_gt = np.zeros_like(b)
        elif nb_ref == 0:
            aa = np.array([1, 1])
            decoded_gt = np.tile(aa, self.ds.shape[1]).reshape((1, self.ds.shape[1], 2))
        elif nb_alt == 2:
            decoder = lambda x: [1, -1] if np.all(x == 2) else [0, 0]
            # np.all() because of b.shape
            decoded_gt = np.apply_along_axis(decoder, axis=1, arr=b)
        elif nb_ref == 2:  # symmetric case with ref allele in only 2 pools: individual is RR or RA
            decoder = lambda x: [0, -1] if np.all(x == 2) else [1, 1]
            decoded_gt = np.apply_along_axis(decoder, axis=1, arr=b)
        else:  # nb_alt > 2 and nb_alt < 8: # nb_alt = 2*n with n > 1
            decoded_gt = np.apply_along_axis(self.multidecoder_gt, axis=1, arr=b)

        return decoded_gt

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def decode_genotypes_gp(self, pooled):
        """
        Recomputes genotypes likelihoods of samples with/without pooling/missing data
        :param samples_gl: samples' true genotypes with phase
        :param dict_gl: likelihoods values to set when encountering missing genotypes
        :return: individual samples genotypes (genotype likelihoods)
        """
        # cdef double[:] scores, rowcolcounts
        # cdef double[:, :] colcross, masks, crosses, unknown, decoded_gl
        # cdef np.ndarray[double, ndim=2] keys

        scores = pooled.flatten()  # np.apply_along_axis(sum, axis=-1, arr=pooled).flatten()
        rowcolcounts = self.rowcolcounter(scores)
        colcross = np.apply_along_axis(np.multiply, 1, self.ds.T, scores)  # (16, 8)
        masks = np.ma.masked_where(self.ds.T, colcross)
        # Returns and sorts genotypes of pools cross section for an individual
        crosses = np.sort(colcross[masks.mask].reshape((self.ds.shape[1], 2)))  # (self.ds.shape[1], 2) = (16, 2) and 2 is the weight of the design
        # this sorts colcross coordinates only (as sorted in the csv table too)
        keys = np.asarray([[*rowcolcounts, *crs] for crs in crosses], dtype=int)
        unknown = np.apply_along_axis(self.multidecoder_gp, axis=1, arr=keys)
        decoded_gl = np.asarray(unknown)

        return decoded_gl

    @staticmethod
    @cython.boundscheck(False)
    @cython.wraparound(False)
    def rowcolcounts(a):
        """
        Count number of pooled RR|RA|AA genotypes over all rows and columns of a pooled matrix
        :param a: score i.e. trinary-encoded true genotypes
        :return: counts of genotypes for the rows and columns
        """
        # cdef double[:] count_rr, count_aa, count_ra, rowcolcounts

        # a should be scores
        count_rr = np.where(a == 0, 1, 0)
        count_aa = np.where(a == 2, 1, 0)
        count_ra = np.where(a == 1, 1, 0)

        rowcolcounts = np.zeros((3*2,), dtype=float)  # counts number of RR, RA, AA rows, same for columns
        rowcolcounts[0] = np.sum(count_rr[:4])
        rowcolcounts[3] = np.sum(count_rr[4:])
        rowcolcounts[1] = np.sum(count_ra[:4])
        rowcolcounts[4] = np.sum(count_ra[4:])
        rowcolcounts[2] = np.sum(count_aa[:4])
        rowcolcounts[5] = np.sum(count_aa[4:])

        return rowcolcounts

    @staticmethod
    @cython.boundscheck(False)
    @cython.wraparound(False)
    def multidecoder_gt(a):
        """
        Decodes pooled scores into individual GT.
        :param a: score
        :return: true genotype with phase
        """

        if np.all(a == 2.0):  # RA * RA
            gt = np.asarray([-1, -1], dtype=float)
        elif np.all(a == 1.0) or np.all(a == 0.0):  # RA * RR or RR * RR
            gt = np.asarray([0, 0], dtype=float)
        else:
            gt = np.asarray([1, 1], dtype=float)
        return gt

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def multidecoder_gp(self, k):

        dkey = get_dummy_key(k)
        gidx = np.digitize(dkey.dot(self.D.transpose()), [len(k)])
        gp = gidx.dot(self.V)
        return gp


# class Decoder(object):
#     """
#         Simulate edeoding step in pooling.
#         Proceed block-wise.
#         This version is based on the use of dot-products of NumPy arrays for
#          representing the lookup table with the adaptive GP values.
#         """
#
#     def __init__(self, design_matrix: np.ndarray, lookup_keys: np.ndarray, lookup_vals: np.ndarray, format: str = 'gt'):
#         self.dm = design_matrix  # matrix for all blocks
#         self.ds1 = Design()  # single block
#         self.dm1 = self.ds1.matrix
#         self.fmt = format.upper()
#         assert (self.fmt == 'GT' or self.fmt == 'GP'), 'Pooling to other formats than GT or GP not implemented'
#         self.D, self.V = lookup_keys, lookup_vals
#         #TODO: if GP and None lookup -> decode into 0.33, 0.33, 0.33
#
#
#     @property
#     def n_blocks(self):
#         return self.dm.shape[1] // self.dm1.shape[1]
#
#     def decode_genotypes_gt(self, pooled):
#         """
#         Recomputes true genotypes of samples with/without pooling/missing data
#         :param pooled: sum of alleles of pooled true genotypes (unpohased)
#         :return: individual samples genotypes (true genotype unphased)
#         """
#         where_alt = np.where(pooled >= 1, 1, 0)
#         where_ref = np.where(pooled <= 1, 1, 0)
#
#         count_alt = where_alt.reshape((1, self.n_blocks, self.dm1.shape[0])).sum(axis=-1)  # 1 count per block
#         count_ref = where_ref.reshape((1, self.n_blocks, self.dm1.shape[0])).sum(axis=-1)
#
#         nb_alt = np.repeat(count_alt, self.dm1.shape[1])  # 1 count per sample
#         nb_ref = np.repeat(count_ref, self.dm1.shape[1])
#
#         encoded = np.dot(pooled, self.dm).reshape(1, self.dm.shape[1], 1)
#         scores = np.full((5, self.dm.shape[1]), -1, dtype=int)
#         # 5 rows = 1 row for nb_alt, 1 row for nb_ref, 1 row for encoded, 2 rows for the decoded genotype
#         scores[0] = nb_alt
#         scores[1] = nb_ref
#         scores[2] = encoded.squeeze()
#         decoded_gt = np.apply_along_axis(self.multidecoder_gt, axis=0, arr=scores)[3:]
#
#         return decoded_gt.T
#
#     @staticmethod
#     def multidecoder_gt(a):
#         """
#         Decodes pooled scores into individual GT.
#         :param a: score as quintet
#         :return:
#         """
#         if a[0] == 0:  # all RR
#             a[3] = 0
#             a[4] = 0
#         elif a[1] == 0:  # all AA
#             a[3] = 1
#             a[4] = 1
#         elif a[0] == 2:   # all RR but 1 AA or RA
#             if a[2] == 2:
#                 a[3] = 1
#             else:
#                 a[3] = 0
#                 a[4] = 0
#         elif a[1] == 2:   # symmetric case: all AA but 1 RR or RA
#             if a[2] == 2:
#                 a[3] = 0
#             else:
#                 a[3] = 1
#                 a[4] = 1
#         else:  # mix of RR, AA, RA/AR
#             if a[2] > 2:
#                 a[3] = 1
#                 a[4] = 1
#             elif a[2] < 2:
#                 a[3] = 0
#                 a[4] = 0
#         return a
#
#     def decode_genotypes_gp(self, pooled):
#         """
#
#         """
#         # Reshape for NumPy style np.apply_along_axis
#         scores = pooled.reshape((self.n_blocks, self.dm1.shape[0]))
#         rowcolcounts = np.apply_along_axis(self.rowcolcounts, axis=-1, arr=scores)
#         # Gets scores of the pair of pools intersecting at each sample
#         colcross = np.apply_along_axis(np.multiply, 1, self.dm.T, pooled).squeeze()
#         # Masks non-zeros i.e. pair of pools scores
#         masks = np.ma.masked_where(self.dm.T, colcross)
#         # Returns and sorts genotypes of pools cross section for an individual
#         crosses = np.sort(colcross[masks.mask].reshape((self.dm.shape[1], 2)))  # 2 is the weight of the design
#         # this sorts colcross coordinates only (as sorted in the csv table too)
#         poolscounts = np.tile(rowcolcounts,
#                               self.dm1.shape[1]).reshape(self.dm.shape[1], rowcolcounts.shape[1])
#         keys = np.concatenate([poolscounts, crosses], axis=1)
#         decoded_gp = np.apply_along_axis(self.multidecoder_gp, axis=-1, arr=keys)
#
#         return decoded_gp
#
#     @staticmethod
#     def rowcolcounts(a):
#         """
#         Count number of pooled RR|RA|AA genotypes over all rows and columns of a pooled matrix
#         :param a: score i.e. trinary-encoded true genotypes
#         :return: counts of genotypes for the rows and columns
#         """
#         # a should be scores
#         count_rr = np.where(a == 0, 1, 0)
#         count_aa = np.where(a == 2, 1, 0)
#         count_ra = np.where(a == 1, 1, 0)
#
#         rowcolcounts = np.zeros((3*2,), dtype=int)  # counts number of RR, RA, AA rows, same for columns
#         rowcolcounts[0] = np.sum(count_rr[:4])
#         rowcolcounts[3] = np.sum(count_rr[4:])
#         rowcolcounts[1] = np.sum(count_ra[:4])
#         rowcolcounts[4] = np.sum(count_ra[4:])
#         rowcolcounts[2] = np.sum(count_aa[:4])
#         rowcolcounts[5] = np.sum(count_aa[4:])
#
#         return rowcolcounts
#
#     def multidecoder_gp(self, k):
#         dkey = get_dummy_key(k)
#         gidx = np.digitize(dkey.dot(self.D.transpose()), [len(k)])
#         gp = gidx.dot(self.V)
#         return gp


def load_lookup_table(path):
    """
    Provides adaptive GL values as a DataFrame
    """
    df = pd.read_csv(os.path.join(path, 'adaptive_gls.csv'),
                     header=None,
                     names=['rowsrr', 'rowsra', 'rowsaa', 'colsrr', 'colsra', 'colsaa',
                            'n', 'm',
                            'rr', 'ra', 'aa']
                     )
    return df


def get_lookup_arrays(path):
    """
    Converts a lookup table to:
    * D: a categorical array encoding key as categorical ("dummy" array)
    * V: a value array where the position of any value matches its key position in D
    """
    T = load_lookup_table(path)
    T.drop_duplicates(inplace=True)
    T.reset_index(drop=True, inplace=True)

    dumlist = []
    V = T[T.columns[-3:]]
    for col in T.columns[:-3]:
        dumlist.append(pd.get_dummies(T[col], prefix=col))
    D = dumlist.pop(0)
    for df in dumlist:
        D = D.join(df)

    D = D.values.astype(float)
    V = V.values.astype(float)
    return D, V


def get_dummy_key(k) :
    """
    Converts the key array to a categorical "dummy" array
    """
    ds = Design()
    rg = ds.pools_size + 1
    strides = np.asarray([rg, rg, rg, rg, rg, rg, 3, 3])  # 3 = genotypes RR, RA/AR, AA
    dumkey = np.zeros((strides.sum(),), dtype=int)
    idx = 0
    for i in range(len(k)):
        dumkey[idx + int(k[i])] = 1
        idx = idx + strides[i]
    return dumkey


def load_lookup_dict(path):
    df = pd.read_csv(os.path.join(path, 'adaptive_gls.csv'),
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
    return df2dict


def blocks_decoder(nB, v, step, lookup_keys, lookup_vals, dec_fmt: str):
    """
    Decodes a single block from a NORB pooling design
    """
    ds = Design()
    dm = ds.matrix
    decoder = SingleBlockDecoder(dm, lookup_keys, lookup_vals, dec_fmt)
    res = []
    cdef int b
    for b in range(nB):
        p = v.squeeze()[step*b:step*b + step]
        if decoder.fmt == 'GP':
            q = decoder.decode_genotypes_gp(p)
        elif decoder.fmt == 'GT':
            q = decoder.decode_genotypes_gt(p)
        res.append(q)
    return np.asarray(res).squeeze()


def single_block_decoder(p, lookup_keys, lookup_vals, dec_fmt):
    ds = Design()
    dm = ds.matrix
    decoder = SingleBlockDecoder(dm, lookup_keys, lookup_vals, dec_fmt)
    if decoder.fmt == 'GP':
        q = decoder.decode_genotypes_gp(p.astype(float))
    elif decoder.fmt == 'GT':
        q = decoder.decode_genotypes_gt(p.astype(float))
    return q

