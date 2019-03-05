import os
import argparse
import subprocess
import warnings
import numpy as np
import time

from scipy.stats import bernoulli as bn
from cyvcf2 import VCF, Writer
from operator import *

from scripts.poolSNPs import parameters as prm
from persotools.debugging import *
from persotools.files import *

warnings.simplefilter('ignore')

global source, KIND, MSS, POOL, SAMPLES, CHK_SZ, WD

### README
#TODO: remove miss.pool
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


def split_pools(idv_id, pool_size, seed=None):
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
    if seed != None:
        np.random.seed(seed)
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
        cls.samples = list()
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

    def set_subset(self, subset): #subsamp from split_pools
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

    def get_subset(self):
        ids = self.flatten()#.reshape(1, self.size)
        return ids

    def pools_list(self):
        design = self.design_matrix()
        if np.where(self == '', False, True).all():
            pools_ = []
            for i in range(design.shape[0]):
                cache = (design[i,:].reshape(self.shape) == False)
                pool = np.ma.masked_array(self, mask=cache)
                pools_.append(pool.compressed())
            return pools_
            # return just for being able to print list if wished

    def get_line(self, samples, variant, kind):
        self.kind = kind
        self.variant = variant.genotypes
        self.samples = samples

    def discretize_genotypes(self):
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
        subset = self.get_subset()
        idx = np.argwhere(np.isin(self.samples, subset))
        self.call = np.asarray(self.variant)[idx]
        return self.call

    def pool_genotypes(self):
        """
        Computes genotypes of the different pools.
        :return: array of {0,1} GLs for RR|RA|AA for each pool
        """
        call = self.get_call().reshape((1, self.size, 3))
        scores = np.apply_along_axis(sum, axis=-1, arr=call[:, :, :-1])

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
            pooler = lambda x: ([0, 0, 0] if np.all(x == 0)
                                else ([1, 1, 0] if np.all(x == self.pools_nb)
                                      else [1, 0, 0]))
            p = np.apply_along_axis(pooler, axis=-1, arr=pooled_gt)

        return p  # list of gt for the 8 pools from design matrix

    def decode_genotypes(self, samples_gt, drop=False):
        """
        Recoomputes genotypes of samples with/without pooling/missing data
        :param pooled_samples: Variant.genotypes
        :return:
        """
        pooled = self.pool_genotypes()
        scores = np.apply_along_axis(sum, axis=-1, arr=pooled[:, :, :-1])
        p = np.argwhere(np.isin(self.samples, self.subset))

        nb_alt = np.sum(
            np.apply_along_axis(
                lambda x: 1 if 1 in x else 0, axis=0, arr=pooled[:, :, :-1]
            )
        )

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
            elif nb_alt == 2:
                decoder = lambda x: [1, -1, 0] if np.all(x == 2) else [0, 0, 0]
                # np.all() because of b.shape
                decoded_gt = np.apply_along_axis(decoder, axis=-1, arr=b)
            else:  # nb_alt > 2 and nb_alt < 8: # nb_alt = 2*n with n > 1
                decoder = lambda x: ([-1, -1, 0] if np.all(x == 2)
                                     else ([0, 0, 0] if (np.all(x == 1) or np.all(x == 0))
                                           else [1, 1, 0]))
                decoded_gt = np.apply_along_axis(decoder, axis=-1, arr=b)
        np.put_along_axis(samples_gt,
                          np.broadcast_to(p, (self.size, 3)),
                          decoded_gt.squeeze(),
                          axis=0)

        if drop == True:
            #TODO: move to pool_genotypes
            # missing data simulation
            np.put(samples_gt, p, [-1, -1, 0])

        return samples_gt

    def __array_finalize__(self, obj):
        """
        Constructor needed for subclassing NumPy arrays.
        See online documentation.
        :param obj:
        :return:
        """
        if obj is None: return
        self.info = getattr(obj, 'info', None)


#TODO: specific processing for the remaining samples?


### PROCESS SUMMARY

def get_data_chunk():
    pass


def process_file(data, groups, f, fileout):
    """
    Computes and rewrites genotypes of all individuals for all samples.
    :param f: integer, index of the file to process in the list
    """
    print('Missing data: ', MSS[f])
    print('Pooling: ', POOL[f])
    w = Writer(fileout[f], data)
    w.set_threads(4)
    tm = time.time()
    #for n, variant in enumerate(data('20:55167111-55167111')):
    for n, variant in enumerate(data):
        process_line(groups, f, w, variant)
        if n%1000 == 0:
            print('{} variants processed in {:06.2f} sec'.format(n+1, time.time()-tm).ljust(80, '.'))
        # if n+1 == 1000:
        #     break


def process_line(groups, f, w, v, write=True):
    """
    Computes and rewrites genotypes of all individual for one sample.
    :param f: integer, index of the file to process in the list
    :param v: cyvcf2.cyvcf2.Variant object
    :param w: cyvcf2.cyvcf2.Writer object
    :param cnt: integer, counts missing data
    :return: variant object with pooled values for GT/GL
    """
    var = v # copy the variant object to make it readable several times
    pooled_samples = np.asarray(var.genotypes)
    sets = []
    for gp in groups[0]:
        sets.append(SNPsPool().set_subset(gp))
    if POOL[f]:
        for p in sets:
            p.get_line(SAMPLES, var, KIND)
            dlt = random_delete(activate=MSS[f])
            pooled_samples = p.decode_genotypes(pooled_samples, drop=bool(dlt))
    else:
        for p in sets:
            p.get_line(SAMPLES, var, KIND)
            dlt = random_delete(activate=MSS[f])
            idx = np.argwhere(np.isin(SAMPLES, p))
            if dlt:
                np.put(pooled_samples, idx , np.asarray([-1, -1, 0]))

    var.genotypes = pooled_samples.tolist()
    if write:
        w.write_record(var)
    else:
        return var


def init_chunk(chunk=True):
    """

    :param chunk: boolean. If True, a new subfile from 1000
    randomly drawm markers is created i.e. a new set of SNPs.
    :return: list of 16-sized groups of samples i.e. pools
    """
    PATH_IN = 'ALL.chr20.snps.gt.vcf.gz'

    os.chdir(WD)
    raw = VCF(PATH_IN, threads=4)  # VCF iterator object
    splits = split_pools(raw.samples, 16, seed=123)  # list of lists
    with open('ALL.chr20.snps.gt.allID.txt', 'w') as f:
        for s in splits[0]:
            for i in s:
                f.write(i + os.linesep)

    if chunk:
        delete_file('ALL.chr20.snps.gt.chunk{}.vcf.gz'.format(str(CHK_SZ)))
        delete_file('chunk_{}.vcf'.format(str(CHK_SZ)))
        delete_file('TMP.chr20.snps.gt.chunk.vcf')

        print('Extracting header'.ljust(80, '.'))
        subprocess.run('bcftools view -h -Ov -o headers.ALL.chr20.snps.gt.vcf ALL.chr20.snps.gt.vcf.gz',
                       shell=True,
                       cwd=WD)
        print('Sampling markers'.ljust(80, '.'))
        subprocess.run('bcftools view -H ALL.chr20.snps.gt.vcf.gz '
                       + '| sort -R | head -{} > chunk_{}.vcf'.format(str(CHK_SZ), str(CHK_SZ)),
                       shell=True,
                       cwd=WD)
        subprocess.run('cat headers.ALL.chr20.snps.gt.vcf chunk_{}.vcf '.format(str(CHK_SZ))
                       + '> TMP.chr20.snps.gt.chunk.vcf',
                       shell=True,
                       cwd=WD)
        print('BGzipping chunk file'.ljust(80, '.'))
        subprocess.run('bcftools view -Oz -o ALL.chr20.snps.gt.chunk{}.vcf.gz '.format(str(CHK_SZ))
                       + 'TMP.chr20.snps.gt.chunk.vcf',
                       shell=True,
                       cwd=WD)
        delete_file('chunk_{}.vcf'.format(str(CHK_SZ)))
        delete_file('TMP.chr20.snps.gt.chunk.vcf')

        # Eliminate the remaining samples (not involved in any pool)
        subprocess.run(' '.join(['bcftools view -Ov -o ALL.chr20.snps.gt.chunk{}.vcf'.format(str(CHK_SZ)),
                                 '-s',
                                 '^' + ','.join(splits[1]),
                                 'ALL.chr20.snps.gt.chunk{}.vcf.gz'.format(str(CHK_SZ))]),
                       shell=True,
                       cwd=WD)
        subprocess.run('bcftools view -Oz -o ALL.chr20.snps.gt.chunk{}.vcf.gz '.format(str(CHK_SZ))
                       + 'ALL.chr20.snps.gt.chunk{}.vcf'.format(str(CHK_SZ)),
                       shell=True,
                       cwd=WD)
        subprocess.run('bcftools index -f ALL.chr20.snps.gt.chunk{}.vcf.gz'.format(str(CHK_SZ)),
                       shell=True,
                       cwd=WD)
        subprocess.run('bcftools sort -Oz -o ALL.chr20.snps.gt.chunk{}.vcf.gz '.format(str(CHK_SZ))
                       + 'ALL.chr20.snps.gt.chunk{}.vcf.gz'.format(str(CHK_SZ)),
                       shell=True,
                       cwd=WD)
        subprocess.run('bcftools index -f ALL.chr20.snps.gt.chunk{}.vcf.gz'.format(str(CHK_SZ)),
                       shell=True,
                       cwd=WD)

    return splits


def run(splits, iteration):
    os.chdir(WD)

    for it in iteration:
        vcf = VCF(source, threads=4)
        print('Iteration --> ', it)
        start = time.time()
        process_file(vcf, splits, it, PATH_OUT)
        stop = time.time()
        print('Elapsed time for pooling the VCF files: {:06.2f} sec'.format(stop - start))


def write_truncate_vcf(path_in, path_out, trunc):
    w = Writer(path_out, VCF(path_in))
    for i, v in enumerate(VCF(path_in)):
        if i == trunc:
            break
        else:
            w.write_record(v)


if __name__ == '__main__':
    WD = prm.WD
    os.chdir(WD)
    CHK_SZ = prm.CHK_SZ
    PATH_IN = prm.PATH_IN
    PATH_OUT = prm.PATH_OUT
    KIND = prm.KIND
    MSS = prm.MSS
    POOL = prm.POOL
    source = prm.SOURCE

    items = init_chunk(chunk=False)

    data = VCF(source)
    SAMPLES = data.samples

    #run(items, range(2))

    for fi in PATH_OUT:
        trc = prm.SUBCHUNK
        write_truncate_vcf(fi, fi.replace('chunk' + str(CHK_SZ), 'chunk' + str(trc)), trc)
    print('Creating subchunk file'.ljust(80, '.'))
    write_truncate_vcf(source[:-3], source[:-3].replace('chunk' + str(CHK_SZ), 'chunk' + str(trc)), trc)
    # source[:-3] deleting '.gz' in the name
    print('BGzipping subchunk file'.ljust(80, '.'))
    subprocess.run('bcftools view -Oz -o ALL.chr20.snps.gt.chunk{}.vcf.gz '.format(str(trc))
                   + 'ALL.chr20.snps.gt.chunk{}.vcf '.format(str(trc)),
                   shell=True,
                   cwd=WD)
    print('Sorting and indexing subchunk file'.ljust(80, '.'))
    subprocess.run('bcftools sort -Oz -o ALL.chr20.snps.gt.chunk{}.vcf.gz '.format(str(trc))
                   + 'ALL.chr20.snps.gt.chunk{}.vcf.gz'.format(str(trc)),
                   shell=True,
                   cwd=WD)
    subprocess.run('bcftools index -f ALL.chr20.snps.gt.chunk{}.vcf.gz'.format(str(trc)),
                   shell=True,
                   cwd=WD)

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