import os
import argparse
import subprocess
import warnings
import numpy as np
import time
import itertools
import multiprocessing as mp
from scipy.stats import bernoulli as bn
from cyvcf2 import VCF, Writer
from operator import *

warnings.simplefilter('ignore')

global DATA, KIND, MSS, POOLED, SAMPLES, SETS

### README

"""
README: 
Requirements: 
- bcftools

Usage for creating pooled genotypes likelihoods
from a set of samples in a VCF file:
(1) Instanciate a Pool structure:
    pool = Pool()
(2) Update the pool's parameters with file's data:
    pool.get_data(<path_to_file>)
(3) Fill the pooling matrix with data from file:
    pool.fill_in()
(4) Implement the needed class methods
    ex. pool.pools_list()
    
Running from cmd:
$ python3 pool.py /home/camilleclouard/PycharmProjects/1000Genomes/data/tests-beagle/ALL.chr20.snps.gt.vcf.gz gt False /home/camilleclouard/PycharmProjects/1000Genomes/data/tests-beagle/ALL.chr20.pooled.snps.gt.vcf.gz
"""

### TOOLS FUNCTIONS


def delete_file(file_path):
    """

    :param file_path:
    :return:
    """
    if os.path.exists(file_path):
        os.remove(file_path)
    else:
        print("The file does not exists")


def file_size(fname):
    """

    :param fname:
    :return:
    """
    statinfo = os.stat(fname)
    return statinfo.st_size


def random_delete(rate=0.01, activate=False):
    """

    :param arr: variant.genotypes instance
    :param rate:
    :param activate:
    :return:
    """
    if activate:
        flag = bn.rvs(p=rate, size=1)
    else:
        flag=False

    return flag


### CLASS AND METHODS FOR POOLING A SET OF SAMPLES


def split_pools(idv_id, pool_size):
    """

    :param idv_id:
    :param pool_size:
    :return:
    """
    idv_nb = len(idv_id)
    nb_pool, left = divmod(idv_nb, pool_size)
    print('Number of {}-sized pools to create, number of single samples remaining: '.format(pool_size),
          nb_pool, left)
    # Create random pools
    np.random.seed(123)
    pools = np.random.choice(idv_id, size=(nb_pool, pool_size), replace=False)
    left = np.isin(idv_id, pools.flatten())
    singles = np.extract(~left, idv_id)
    pools = pools.tolist()

    return pools, singles


class SNPsPool(np.ndarray):
    """
    Simulates the different steps of a pooling process and
    builds the pooling design.
    """
    def __new__(cls, shape=(4,4), id_len=8, pools_nb=8, pools_size=4):
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
        cls.subset = list()
        cls.samples = dict()
        cls.variant = object
        cls.call = []
        cls.kind = ''
        cls.miss = False
        return np.empty_like(super(SNPsPool, cls).__new__(cls, shape),
                             dtype=id)

    def design_matrix(self, random=False):
        """
        That function is not intended to be called explicitly.
        :param random: bool for dispatching idv randomly in the matrix?
        :return: design matrix
        """
        pools_size = self.pools_size
        design =  np.zeros((self.pools_nb, self.size), dtype=int)
        if random == False:
            for i in range(int(self.pools_nb/self.ndim)):
                j = i * pools_size
                design[i, j:j+pools_size] = [1]*pools_size
            for i in range(int(self.pools_nb/self.ndim), self.pools_nb):
                j = i - pools_size
                design[i,
                       [j+k*pools_size for k in range(pools_size)]] = 1
        return design

    def fill_in(self, subset): #subsamp from split_pools
        """
        Fills the pooling matrix according to the (default) design
        and the input list of samples.
        :param subset: 1D-nparray-like object with the variables IDs
        :return: pooling matrix with samples' names.
        """
        self.subset = subset
        sub = self.subset
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

    def pools_list(self): #after fill_in
        design = self.design_matrix()
        if np.where(self == '', False, True).all():
            pools_ = []
            for i in range(design.shape[0]):
                cache = (design[i,:].reshape(self.shape) == False)
                pool = np.ma.masked_array(self, mask=cache)
                pools_.append(pool.compressed())
            return pools_
            # return just for being able of printing list if wished

    def get_line(self, samples, variant, kind):
        self.kind = kind
        self.variant = variant.genotypes
        self.samples = samples

    def overw_G(self):
        """
        Compute binary genotypes for homozygous individuals,
        the heterozygous ones are set to NaN/unknown.
        :param filename: vcf with subset of idv
        :return: array of {0,1} GLs for RR|RA|AA for each sample
        """
        if self.kind == 'gl':
            pass
            #TODO: implement
            bin_gl = self.variant
            for var in self.variant[:,:-1]:
                lkh = np.vectorize(lambda x: pow(10, x))
                gl = np.around(lkh(var)).astype(float)
                unknown = np.where(np.apply_along_axis(sum, 2, gl) == 0.0)
                for i, j in zip(unknown[0], unknown[1]):
                    gl[i,j,:] = np.full((3,), np.nan)
                return self.variant

    def get_call(self):
        self.call = list() + [self.variant[self.samples[s]] for s in self.subset]
        return np.asarray(self.call)

    def sum_up_G(self, idvs, frompool=False):
        """

        :param idvs: array of samples with genotypes
        :return: pooled 'genotype' over the input array
        """
        idvs = np.asarray(idvs)
        if self.kind == 'gl':
            #TODO: implement
            if [0., 0., 1.] in idvs:
                pooled_G = np.asarray([0, 0., 1.])
            else:
                pooled_G = np.asarray([1., 0., 0.])

        else:
            if np.sum(idvs[:, :-1])//2 == len(idvs):
                pooled_G = np.asarray([1, 1, 0])
            elif np.any(np.isin(idvs[:, :-1], 1)):
                pooled_G = np.asarray([1, -1, 0])
            elif np.any(np.isin(idvs[:, :-1], -1)):
                pooled_G = np.asarray([-1, -1, 0])
            else:
                pooled_G = np.asarray([0, 0, 0])

        return pooled_G # REF, ALT but no phase

    def pools_G(self):
        """
        Computes genotypes of the different pools.
        :return: array of {0,1} GLs for RR|RA|AA for each pool
        """
        design = self.design_matrix()
        call = self.get_call()
        if np.where(self == '', False, True).all():
            pools_G = []
            for i in range(design.shape[0]):
                cache = design[i, :]
                ixgrid = np.ix_([True]*call.shape[0], cache)
                poolval = call[ixgrid]
                pools_G.append(self.sum_up_G(poolval))

        return pools_G # list of gt for the 8 pools from design matrix

    def new_sample_G(self, pooled_samples, drop=False):
        """
        Recoomputes genotypes of samples with/without pooling/missing data
        :param pooled_samples: Variant.genotypes
        :return:
        """
        pools_list = self.pools_list()
        pooled = self.pools_G()
        for s in self.subset: #16
            p_p = np.where(np.isin(np.asarray(pools_list), s))[0]
            # np.put(pooled_samples[s],
            #        [0, 1],
            #        self.sum_up_G(np.asarray([pooled[p] for p in p_p]),
            #                      frompool=True)
            #        )
            pooled_samples[s] = self.sum_up_G(np.asarray([pooled[p] for p in p_p]),
                                 frompool=True)
            if drop == True:
                np.put(pooled_samples[s],
                       [0, 1],
                       [-1, -1]
                       )

        return pooled_samples

    def __array_finalize__(self, obj):
        """
        Constructor needed for subclassing NumPy arrays.
        See online documentation.
        :param obj:
        :return:
        """
        if obj is None: return
        self.info = getattr(obj, 'info', None)


### ARGUMENTS PARSING AND EXECUTION
"""
parser = argparse.ArgumentParser(description='''Computes pooled GL/GT data for each subset.
 Subsets are written as new VCF files''')
parser.add_argument('file_in',
                    action='store',
                    help='VCF file to process')
parser.add_argument('kind',
                    action='store',
                    help='gl/gt')
parser.add_argument('missing',
                    action='store',
                    help='Boolean for simulating missing data (default proportion = 1%)')
parser.add_argument('file_out',
                    action='store',
                    help='VCF file processed')
args = parser.parse_args()
"""
"""
pll = mp.Pool(processes=4)
params = list(zip(files_list, [args.kind]*len(files_list)))
_ = pll.starmap(process_pool, params)
"""
#TODO: specific processing for the remaining samples


### PROCESS SUMMARY

def get_data_chunk():
    chunk = []
    n = 0
    for variant in DATA:
        chunk.append(variant)
        n+=1


def process_file(f):
    """
        Computes and rewrites genotypes of all individuals for all samples.
        :param f: integer, index of the file to process in the list
    """
    print('Missing data: ', MSS[f])
    w = Writer(PATH_OUT[f], DATA)
    tm = time.time()
    cnt = 0
    for n, variant in enumerate(DATA):
        cnt = process_line(f, w, variant, cnt)
        if n%1000 == 0:
            print('{} variants processed in {:06.2f} sec'.format(n+1, time.time()-tm).ljust(80, '.'))
        if n+1 == 1000:
            break
    #DATA.add_to_header('##Number of missing genotypes: {}'.format(cnt))
    #DATA.add_to_header('##Percentage of missing genotypes: {}'.format(cnt*100/(n*len(DATA.samples))))

def process_line(f, w, v, cnt):
    """
    Computes and rewrites genotypes of all individual for one sample.
    :param f: integer, index of the file to process in the list
    :param v: cyvcf2.cyvcf2.Variant object
    :param w: cyvcf2.cyvcf2.Writer object
    :param cnt: integer, counts missing data
    :return: variant object with pooled values for GT/GL
    """
    #tm = time.time()
    var = v # copy the variant object to make it readable several times
    pooled_samples = dict(zip(SAMPLES.keys(), np.array(var.genotypes)))
    if POOLED[f]:
        for p in SETS:
            p.get_line(SAMPLES, var, KIND)
            dlt = random_delete(activate=MSS[f])
            pooled_samples = p.new_sample_G(pooled_samples, drop=bool(dlt))
    else:
        for p in SETS:
            p.get_line(SAMPLES, var, KIND)
            dlt = random_delete(activate=MSS[f])
            if dlt:
                for s in np.array(p).flatten():
                    np.put(pooled_samples[s],
                           [0, 1],
                           [-1, -1]
                           )
    output = sorted([[k, SAMPLES[k], pooled_samples[k]] for k in SAMPLES.keys()], key=itemgetter(1))
    var.genotypes = [out[2] for out in output]
    mask = np.ma.masked_where(np.array(var.genotypes)[0,1].sum() >= 0,
                              np.array(var.genotypes))
    cnt += mask.count()
    w.write_record(var)
    #print('Variant {} processed in {:06.2f} sec'.format(var.ID, time.time() - tm).ljust(80, '.'))
    return cnt


if __name__ == '__main__':
    # np.random.seed(789)  # fixed seed: same output each time, in split func directly

    wd = '/home/camilleclouard/PycharmProjects/1000Genomes/data/tests-beagle'
    os.chdir(wd)
    PATH_IN = 'ALL.chr20.snps.gt.vcf.gz'
    PATH_OUT = ['ALL.chr20.pooled.snps.gt.chunk.vcf',
                'ALL.chr20.pooled.missing.snps.gt.chunk.vcf',
                'ALL.chr20.missing.snps.gt.chunk.vcf']
    raw = VCF(PATH_IN)  # VCF iterator object
    splits = split_pools(raw.samples, 16)  # list of lists

    # iteration = 0
    # iteration = 1
    iteration = 2

    if iteration == -1:
        delete_file('ALL.chr20.snps.gt.chunk.vcf.gz')
        delete_file('chunk_1000.vcf')
        delete_file('TMP.chr20.snps.gt.chunk.vcf')

        print('Extracting header'.ljust(80, '.'))
        subprocess.run('bcftools view -h -Ov -o headers.ALL.chr20.snps.gt.vcf ALL.chr20.snps.gt.vcf.gz',
                       shell=True,
                       cwd=wd)
        print('Sampling markers'.ljust(80, '.'))
        subprocess.run('bcftools view -H ALL.chr20.snps.gt.vcf.gz | sort -R | head -1000 > chunk_1000.vcf',
                       shell=True,
                       cwd=wd)
        subprocess.run('cat headers.ALL.chr20.snps.gt.vcf chunk_1000.vcf > TMP.chr20.snps.gt.chunk.vcf',
                       shell=True,
                       cwd=wd)
        print('BGzipping chunk file'.ljust(80, '.'))
        subprocess.run('bcftools view -Oz -o ALL.chr20.snps.gt.chunk.vcf.gz TMP.chr20.snps.gt.chunk.vcf',
                       shell=True,
                       cwd=wd)
        delete_file('chunk_1000.vcf')
        delete_file('TMP.chr20.snps.gt.chunk.vcf')

        # Eliminate the remaining samples (not involved in any pool)
        subprocess.run(' '.join(['bcftools view -Ov -o ALL.chr20.snps.gt.chunk.vcf',
                                  '-s',
                                  '^' + ','.join(splits[1]),
                                  'ALL.chr20.snps.gt.chunk.vcf.gz']),
                       shell=True,
                       cwd='/home/camilleclouard/PycharmProjects/1000Genomes/data/tests-beagle')
        subprocess.run('bcftools view -Oz -o ALL.chr20.snps.gt.chunk.vcf.gz ALL.chr20.snps.gt.chunk.vcf',
                       shell=True,
                       cwd='/home/camilleclouard/PycharmProjects/1000Genomes/data/tests-beagle')
        subprocess.run('bcftools index -f ALL.chr20.snps.gt.chunk.vcf.gz',
                       shell=True,
                       cwd='/home/camilleclouard/PycharmProjects/1000Genomes/data/tests-beagle')
        subprocess.run('bcftools sort -Oz -o ALL.chr20.snps.gt.chunk.vcf.gz ALL.chr20.snps.gt.chunk.vcf.gz',
                       shell=True,
                       cwd='/home/camilleclouard/PycharmProjects/1000Genomes/data/tests-beagle')
        subprocess.run('bcftools index -f ALL.chr20.snps.gt.chunk.vcf.gz',
                       shell=True,
                       cwd='/home/camilleclouard/PycharmProjects/1000Genomes/data/tests-beagle')

    DATA = VCF('ALL.chr20.snps.gt.chunk.vcf.gz')
    KIND = 'gt'
    MSS = [False, True, True]
    POOLED = [True, True, False]
    SAMPLES = dict(zip(DATA.samples, range(len(DATA.samples))))  # dict
    oo = SNPsPool()
    SETS = []
    for sp in splits[0]:
        SETS.append(SNPsPool().fill_in(sp))

    design = SNPsPool().design_matrix()

    start = time.time()
    # multiprocessing not usable with the VCF object structure!
    # Smthg goes wrong in the source file when using mp...

    # Run one by one!
    #TODO: cmd version
    process_file(iteration)

    stop = time.time()
    print('Elapsed time for pooling the VCF files: {:06.2f} sec'.format(stop-start))

    ### Mtx test
    # M = np.array([[[1,0,1], [0,0,1], [0,0,1], [0,0,1]],
    #               [[0,0,1], [0,0,1], [0,0,1], [0,0,1]],
    #               [[0,0,1], [0,0,1], [0,0,1], [0,0,1]],
    #               [[0,0,1], [0,0,1], [0,0,1], [0,0,1]]])
    #
    # oo = SNPsPool()
    # p8 = np.zeros((1, 8, 3), dtype=int)
    # for i in range(4):
    #     p = oo.sum_up_G(M[i, :])
    #     p8[:, i, :] = p
    #     print('P{0} --> {1}'.format(i+1, p))
    #     p = oo.sum_up_G(M[:, i])
    #     print('P{0} --> {1}'.format(i+5, p))
    #     p8[:, i+4, :] = p
    #
    # for j in range(4):
    #     for k in range(4,8):
    #         s = oo.sum_up_G(np.concatenate((p8[:, j, :], p8[:, k, :])))
    #         print('Sample I{}{} --> {}'.format(j+1, k+1, s))
    # # pb with all/any logical tests on arrays