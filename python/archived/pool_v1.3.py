import os
import subprocess
import warnings
import numpy as np
import time
from scipy.stats import bernoulli as bn
from cyvcf2 import VCF, Writer
from operator import *

warnings.simplefilter('ignore')

global source, KIND, MSS, POOLED, SAMPLES, chk_sz, wd

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
        numpy.random.seed(seed)
    pools = numpy.random.choice(idv_id, size=(nb_pool, pool_size), replace=False)
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
                pool = numpy.ma.masked_array(self, mask=cache)
                pools_.append(pool.compressed())
            return pools_
            # return just for being able to print list if wished

    def get_line(self, samples, variant):
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
        ids = self.get_subset()
        self.call = [self.variant[self.samples[i]] for i in ids]
        return np.asarray(self.call)

    def pooling_rules(self, idvs):
        """

        :param idvs: array of samples with genotypes
        :return: pooled 'genotype' over the input array
        """
        idvs = np.asarray(idvs)
        if np.sum(idvs[:, :-1])/2 == len(idvs):
            pooled_G = np.asarray([1, 1, 0])
        elif np.any(np.isin(idvs[:, :-1], 1)):
            pooled_G = np.asarray([1, 0, 0])
        elif np.all(np.isin(idvs[:, :-1], -1)):
            pooled_G = np.asarray([-1, -1, 0])
        else:
            pooled_G = np.asarray([0, 0, 0])

        return pooled_G # REF, ALT but no phase

    def decoding_rules(self, pattern, nb_alt):
        # 2 pools involved for each idv
        # --> pattern has 2 elements, 2*2 binary alleles
        pp = pattern[:, :, :-1].squeeze()

        if np.isin(pp, -1).any():
            gt = np.asarray([-1, -1, 0])
        elif np.all(np.isin(pp, 0)):
            gt = np.asarray([0, 0, 0])

        elif nb_alt == 0:
            gt = np.asarray([0, 0, 0])

        elif nb_alt == 1:
            print('Impossible combination. Error occured in pooling.')

        elif nb_alt == 2:
            # 1 idv wih ALT (no matter how many)
            if np.sum(pp) == 2:
                # 1/0 x 1/0
                gt = np.asarray([1, -1, 0])
            elif np.sum(pp) == 1:
                # 0/0 x 1/0
                gt = np.asarray([0, 0, 0])
            else:
                # 1/1 x 1/0
                gt = np.asarray([1, 1, 0])


        else:
            if np.sum(pp) == 0:
                # 0/0 x 0/0
                gt = np.asarray([0, 0, 0])
            elif np.sum(pp) == 1:
                # 0/0 x 1/0
                gt = np.asarray([0, 0, 0])
            elif np.sum(pp) == 2:
                # 1/0 x 1/0
                gt = np.asarray([-1, -1, 0])
            elif np.sum(pp) == 3:
                # 1/1 x 1/0
                gt = np.asarray([1, 1, 0])
            else: # sum = 4
                gt = np.asarray([1, 1, 0])
        return gt


    def pool_genotypes(self):
        """
        Computes genotypes of the different pools.
        :return: array of {0,1} GLs for RR|RA|AA for each pool
        """
        d = self.design_matrix()
        call = self.get_call()
        pools = np.zeros((1, self.pools_nb, 3), dtype=int)
        for i in range(d.shape[0]):
            cache = d[i, :]
            poolval = call[np.argwhere(cache == 1)].squeeze()
            pools[:, i, :] = self.pooling_rules(poolval)
        return pools

    def decode_genotypes(self, pooled_samples, drop=False):
        """
        Recoomputes genotypes of samples with/without pooling/missing data
        :param pooled_samples: Variant.genotypes
        :return:
        """
        pools = self.pools_list()
        pooled = self.pool_genotypes()
        nb_alt = np.sum(
            np.apply_along_axis(
                lambda x: 1 in x,
                -1,
                pooled[:, :, :-1]
            )
        )
        for s in self.get_subset(): #16
            p = np.where(np.isin(np.asarray(pools), s))[0] # pools pair involved
            out = self.decoding_rules(pooled[:, p], nb_alt)
            pooled_samples[s] = out
            if drop == True:
                # missing data simulation
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


#TODO: specific processing for the remaining samples


### PROCESS SUMMARY

def get_data_chunk():
    vcf = VCF(source)
    chunk = []
    n = 0
    for variant in vcf:
        chunk.append(variant)
        n+=1


def process_file(data, groups, f, fileout):
    """
        Computes and rewrites genotypes of all individuals for all samples.
        :param f: integer, index of the file to process in the list
    """
    print('Missing data: ', MSS[f])
    w = Writer(fileout[f], data)
    tm = time.time()
    #for n, variant in enumerate(data('20:55167111-55167111')):
    for n, variant in enumerate(data):
        process_line(groups, f, w, variant)
        if n%1000 == 0:
            print('{} variants processed in {:06.2f} sec'.format(n+1, time.time()-tm).ljust(80, '.'))
        # if n+1 == 1000:
        #     break


def process_line(groups, f, w, v):
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
    sets = []
    for gp in groups[0]:
        sets.append(SNPsPool().set_subset(gp))
    if POOLED[f]:
        for p in sets:
            p.set_line_values(SAMPLES, var, KIND)
            dlt = random_delete(activate=MSS[f])
            pooled_samples = p.decode_genotypes(pooled_samples, drop=bool(dlt))
    else:
        for p in sets:
            p.set_line_values(SAMPLES, var, KIND)
            dlt = random_delete(activate=MSS[f])
            if dlt:
                for s in np.array(p).flatten():
                    np.put(pooled_samples[s],
                           [0, 1],
                           [-1, -1]
                           )
    output = sorted([[k, SAMPLES[k], pooled_samples[k]] for k in SAMPLES.keys()], key=itemgetter(1))
    var.genotypes = [out[2] for out in output]
    w.write_record(var)

def init_chunk(chunk=True):
    """

    :param chunk: boolean. If True, a new subfile from 1000
    randomly drawm markers is created i.e. a new set of SNPs.
    :return: list of 16-sized groups of samples i.e. pools
    """
    wd = '/home/camilleclouard/PycharmProjects/1000Genomes/data/tests-beagle'
    chk_sz = 1000
    PATH_IN = 'ALL.chr20.snps.gt.vcf.gz'

    os.chdir(wd)
    raw = VCF(PATH_IN)  # VCF iterator object
    splits = split_pools(raw.samples, 16, seed=456)  # list of lists

    if chunk:
        delete_file('ALL.chr20.snps.gt.chunk{}.vcf.gz'.format(str(chk_sz)))
        delete_file('chunk_{}.vcf'.format(str(chk_sz)))
        delete_file('TMP.chr20.snps.gt.chunk.vcf')

        print('Extracting header'.ljust(80, '.'))
        subprocess.run('bcftools view -h -Ov -o headers.ALL.chr20.snps.gt.vcf ALL.chr20.snps.gt.vcf.gz',
                       shell=True,
                       cwd=wd)
        print('Sampling markers'.ljust(80, '.'))
        subprocess.run('bcftools view -H ALL.chr20.snps.gt.vcf.gz '
                       + '| sort -R | head -{} > chunk_{}.vcf'.format(str(chk_sz), str(chk_sz)),
                       shell=True,
                       cwd=wd)
        subprocess.run('cat headers.ALL.chr20.snps.gt.vcf chunk_{}.vcf '.format(str(chk_sz))
                       + '> TMP.chr20.snps.gt.chunk.vcf',
                       shell=True,
                       cwd=wd)
        print('BGzipping chunk file'.ljust(80, '.'))
        subprocess.run('bcftools view -Oz -o ALL.chr20.snps.gt.chunk{}.vcf.gz '.format(str(chk_sz))
                       + 'TMP.chr20.snps.gt.chunk.vcf',
                       shell=True,
                       cwd=wd)
        delete_file('chunk_{}.vcf'.format(str(chk_sz)))
        delete_file('TMP.chr20.snps.gt.chunk.vcf')

        # Eliminate the remaining samples (not involved in any pool)
        subprocess.run(' '.join(['bcftools view -Ov -o ALL.chr20.snps.gt.chunk{}.vcf'.format(str(chk_sz)),
                                 '-s',
                                 '^' + ','.join(splits[1]),
                                 'ALL.chr20.snps.gt.chunk{}.vcf.gz'.format(str(chk_sz))]),
                       shell=True,
                       cwd='/home/camilleclouard/PycharmProjects/1000Genomes/data/tests-beagle')
        subprocess.run('bcftools view -Oz -o ALL.chr20.snps.gt.chunk{}.vcf.gz '.format(str(chk_sz))
                       + 'ALL.chr20.snps.gt.chunk{}.vcf'.format(str(chk_sz)),
                       shell=True,
                       cwd='/home/camilleclouard/PycharmProjects/1000Genomes/data/tests-beagle')
        subprocess.run('bcftools index -f ALL.chr20.snps.gt.chunk{}.vcf.gz'.format(str(chk_sz)),
                       shell=True,
                       cwd='/home/camilleclouard/PycharmProjects/1000Genomes/data/tests-beagle')
        subprocess.run('bcftools sort -Oz -o ALL.chr20.snps.gt.chunk{}.vcf.gz '.format(str(chk_sz))
                       + 'ALL.chr20.snps.gt.chunk{}.vcf.gz'.format(str(chk_sz)),
                       shell=True,
                       cwd='/home/camilleclouard/PycharmProjects/1000Genomes/data/tests-beagle')
        subprocess.run('bcftools index -f ALL.chr20.snps.gt.chunk{}.vcf.gz'.format(str(chk_sz)),
                       shell=True,
                       cwd='/home/camilleclouard/PycharmProjects/1000Genomes/data/tests-beagle')

    return splits


def run(splits, iteration):
    os.chdir(wd)

    for it in iteration:
        vcf = VCF(source)
        print('Iteration --> ', it)
        start = time.time()
        process_file(vcf, splits, it, PATH_OUT)
        stop = time.time()
        print('Elapsed time for pooling the VCF files: {:06.2f} sec'.format(stop - start))
        del vcf


if __name__ == '__main__':
    #TODO: encode "hard" pieces in parameters.py
    wd = '/home/camilleclouard/PycharmProjects/1000Genomes/data/tests-beagle'
    os.chdir(wd)
    chk_sz = 1000
    PATH_IN = 'ALL.chr20.snps.gt.vcf.gz'
    PATH_OUT = ['ALL.chr20.pooled.snps.gt.chunk{}.vcf'.format(str(chk_sz)),
                'ALL.chr20.pooled.missing.snps.gt.chunk{}.vcf'.format(str(chk_sz)),
                'ALL.chr20.missing.snps.gt.chunk{}.vcf'.format(str(chk_sz))]
    KIND = 'gt'
    MSS = [False, True, True]
    POOLED = [True, True, False]
    source = 'ALL.chr20.snps.gt.chunk{}.vcf.gz'.format(str(chk_sz))
    data = VCF(source)
    SAMPLES = dict(zip(data.samples, range(len(data.samples))))


    items = init_chunk(chunk = True)
    run(items, range(3))


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