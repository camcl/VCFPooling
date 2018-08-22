import os
import argparse
import subprocess
import warnings
import allel
import numpy as np
import time
import gzip
import gc
import multiprocessing as mp
from scipy.stats import bernoulli as bn
from cyvcf2 import VCF, Writer

warnings.simplefilter('ignore')

### README

"""
README: 
Requirements: 
- AWK/GAWK package has to be installed on the OS
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
"""

### TOOLS FUNCTIONS

def delete_file(file_path):
    if os.path.exists(file_path):
        os.remove(file_path)
    else:
        print("The file does not exists")

def _zip_vcf(file_path):
    """
    Compress a file as a .gz and remove the uncompressed file.
    :param file_path: path to the uncompressed file
    :return:
    """
    inF = open(file_path, 'rb')
    s = inF.read()
    inF.close()

    outF = gzip.GzipFile(file_path + ".gz", 'wb')
    outF.write(s)
    outF.close()

    delete_file(file_path)

def file_size(fname):
    statinfo = os.stat(fname)
    return statinfo.st_size

def random_delete(arr, rate):
    np.random.seed(123)
    flag = bn.rvs(p=rate, size=(arr.shape[0], arr.shape[1])).reshape(arr.shape[0], arr.shape[1])
    missing = -np.ones_like(arr[:, :, 0], dtype=int)
    arr[:, :, 0][flag == 1] = missing[flag == 1]
    arr[:, :, 1][flag == 1] = missing[flag == 1]
    return arr

def process_pool(file_path, kind):
    start = time.time()
    set = SNPsPool()
    set.get_data(file_path, kind)
    set.fill_in()
    set.samples_to_file()
    set.split_vcf()
    set.replace_calls()
    set.reassemble_vcf()
    stop = time.time()
    print('Processing {:d} samples in one batch takes: {:06.2f} sec\r\n'.format(set.size, stop - start))

### SHELL COMMANDS

VCF_HEADER = """
bcftools view -h {0} > {1}
""" # {0} contains vcf file's name to process, and {1} output file with headers

VCF_BODY = """
bcftools view -H {0} > {1}
""" # {0} contains vcf file's name to process, and {1} output file with only calls

### CLASS AND METHODS FOR POOLING A SET OF SAMPLES

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
        cls.path = ''
        cls.samples = list()
        cls.calls = np.zeros((shape[0]*shape[1],), dtype=np.int_)
        cls.data = dict()
        cls.kind = ''
        cls.miss = True
        return np.empty_like(super(SNPsPool, cls).__new__(cls, shape),
                             dtype=id)

    def design_matrix(self, random=False):
        """
        That function is not intended to be called explicitly.
        :param random: bool for dispatching idv randomly in the matrix?
        :return: design matrix
        """
        pools_nb = self.pools_nb
        pools_size = self.pools_size
        design =  np.zeros((pools_nb, self.size), dtype=np.int_)
        if random == False:
            for i in range(int(pools_nb/self.ndim)):
                j = i * pools_size
                design[i, j:j+pools_size] = [1]*pools_size
            for i in range(int(pools_nb/self.ndim), pools_nb):
                j = i - pools_size
                design[i,
                       [j+k*pools_size for k in range(pools_size)]] = 1
        #TODO: case of a random pooling?
        return design

    def get_data(self, path, kind):
        """
        Extract data from the subset VCF.
        :param path: path to the VCF file
        :return: dictionary with lists of samples and their genotypes (GL).
        """
        self.path = path
        self.data = allel.read_vcf(self.path, fields=['samples', 'calldata/' + kind.upper()])
        self.samples = self.data['samples']
        self.calls = self.data['calldata/' + kind.upper()]
        self.kind = kind

    def fill_in(self):
        """
        Fills the pooling matrix according to the (default) design
        and the input list of samples.
        :param subset: 1D-nparray-like object with the variables IDs
        :return: pooling matrix with samples' names.
        """
        subset = list(self.samples)
        try:
            for i in range(self.shape[0]):
                self[i, :] = subset[:self.shape[1]]
                subset = subset[self.shape[1]:]
        except Exception as exc:
            if len(subset) > self.size:
                raise ValueError('The input you gave is too long') from exc
            if len(subset) < self.size:
                raise ValueError('The input you gave is too short') from exc
            if type(subset) != np.ndarray and type(subset) != list:
                raise TypeError('The input is not a 1D-array-like') from exc
            if len(subset) > 0 and type(subset[0]) != str:
                raise TypeError('The input does not contain str-type elements') from exc

    def pools_list(self):
        self.design  = self.design_matrix()
        if np.where(self == '', False, True).all() == True:
            pools_ = []
            #TODO: what if not the default matrix?
            for i in range(self.design.shape[0]):
                cache = (self.design[i,:].reshape(self.shape) == False)
                pool = np.ma.masked_array(self, mask=cache)
                pools_.append(pool.compressed())
            return pools_
            # return just for being able of printing list if wished

    def overw_G(self):
        """
        Compute binary genotypes for homozygous individuals,
        the heterozygous ones are set to NaN/unknown.
        :param filename: vcf with subset of idv
        :return: array of {0,1} GLs for RR|RA|AA for each sample
        """
        if self.kind == 'gl':
            lkh = np.vectorize(lambda x: pow(10, x))
            gl = np.around(lkh(self.calls)).astype(float)
            unknown = np.where(np.apply_along_axis(sum, 2, gl) == 0.0)
            for i, j in zip(unknown[0], unknown[1]):
                gl[i,j,:] = np.full((3,), np.nan)
            return gl

        else:
            return self.calls

    def sim_miss(self):
        """

        :return:
        """
        print('Simulate missing data'.ljust(80, '.'))
        gs = self.overw_G()
        if self.miss == True:
            return random_delete(gs, 0.01)
        else:
            return gs


    def sum_up_G(self, idvs):
        """

        :param idvs: array of samples with genotypes
        :return: pooled 'genotype' over the input array
        """
        idvs = np.asarray(idvs)
        pooled_G = np.zeros((idvs.shape[0], 1, idvs.shape[2]))
        for snp in range(idvs.shape[0]):
            if [np.nan]*idvs.shape[2] in idvs[snp, :, :].tolist():
                pooled_G[[snp], :, :] = np.asarray([np.nan]*idvs.shape[2])
                if self.kind == 'gl':
                    if [0., 0., 1.] in idvs[snp, :, :].tolist():
                        pooled_G[[snp], :, :] = np.asarray([0, 0., 1.])
                    else:
                        pooled_G[[snp], :, :] = np.asarray([1., 0., 0.])

                else:
                    if [1, 1] in idvs[snp, :, :].tolist():
                        pooled_G[[snp], :, :] = np.asarray([1, 1], dtype=int)
                    elif [0, 1] in idvs[snp, :, :].tolist()\
                            and [1, 0] in idvs[snp, :, :].tolist():
                        pooled_G[[snp], :, :] = np.asarray([1, 1], dtype=int)
                    elif [0, 1] in idvs[snp, :, :].tolist() \
                            and [1, 0] not in idvs[snp, :, :].tolist() \
                            and [-1, -1] not in idvs[snp, :, :].tolist():
                        pooled_G[[snp], :, :] = np.asarray([0, 1], dtype=int)
                    elif [0, 1] in idvs[snp, :, :].tolist() \
                            and [1, 0] not in idvs[snp, :, :].tolist() \
                            and [-1, -1] in idvs[snp, :, :].tolist():
                        pooled_G[[snp], :, :] = np.asarray([-1, 1], dtype=int)
                    elif [0, 1] not in idvs[snp, :, :].tolist() \
                            and [1, 0] in idvs[snp, :, :].tolist() \
                            and [-1, -1] in idvs[snp, :, :].tolist():
                        pooled_G[[snp], :, :] = np.asarray([1, -1], dtype=int)
                    elif [0, 1] not in idvs[snp, :, :].tolist() \
                            and [1, 0] in idvs[snp, :, :].tolist() \
                            and [-1, -1] not in idvs[snp, :, :].tolist():
                        pooled_G[[snp], :, :] = np.asarray([1, 0], dtype=int)
                    elif [0, 1] not in idvs[snp, :, :].tolist() \
                            and [1, 0] not in idvs[snp, :, :].tolist() \
                            and [-1, -1] in idvs[snp, :, :].tolist():
                        pooled_G[[snp], :, :] = np.asarray([-1, -1], dtype=int)
                    else:
                        pooled_G[[snp], :, :] = np.asarray([0, 0], dtype=int)

        return pooled_G

    def pools_G(self):
        """
        Computes genotypes of the different pools.
        :return: array of {0,1} GLs for RR|RA|AA for each pool
        """
        gg = self.sim_miss()
        self.design = self.design_matrix()
        if np.where(self == '', False, True).all() == True:
            gg_ = []
            # TODO: what if not the default matrix?
            for i in range(self.design.shape[0]):
                cache = self.design[i, :]
                ixgrid = np.ix_([True]*gg.shape[0], cache, [True]*gg.shape[2])
                gg_.append(gg[ixgrid])
        pools_G = []
        for p in range(len(self.pools_list())):
            pools_G.append(self.sum_up_G(gg_[p]))
        return pools_G

    def new_sample_G(self):
        """
        Simulates GLs from pooling samples together.
        :return: array of {0,1} GLs for RR|RA|AA for each sample
        """
        pooled_samples = np.zeros_like(self.calls)
        pools_list = self.pools_list()
        pooled = self.pools_G()
        for s in range(pooled_samples.shape[1]): #16
            p_p = np.where(np.isin(np.asarray(pools_list), self.samples[s]))[0]
            pooled_samples[:, [s], :] = self.sum_up_G(np.concatenate([pooled[p] for p in p_p],
                                                                     axis=1))
        gc.enable()
        del pooled
        gc.collect()
        return pooled_samples

    def samples_to_file(self):
        """
        Writes old GL values/new ones pooled as mapping columns
        into a tab-delimited file.
        one mapping file per sample.
        :return:
        """
        new_data = self.new_sample_G()
        print('Start writing pooled GL/GT for subset {}'.format(self.path[-10:-7]).ljust(80, '.'))
        with open('new_data_' + self.kind + '_' + self.path[-10:-7] + '.vcf', 'wb') as f:
            #mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_WRITE)
            mm = f
            for row in zip(*[new_data[:, i, :] for i in range(len(self.samples))]):
                if self.kind == 'gl':
                    row = [np.array2string(fld, separator=',')[1:-1].replace(' ', '') for fld in row]
                else:
                    row = [np.array2string(fld, separator='|')[1:-1].replace(' ', '') for fld in row]
                row = '\t'.join(row)
                mm.write((row + '\n').encode())
            #mm.flush()
            #mm.close()

    def split_vcf(self):
        """
        Split the VCF input in headers and body (0on ly calls) sub-vcf.
        Keep on ly the mandatory fixed fields values in the body.
        Preprocessing for replacing the calls' data.
        :return:
        """
        subprocess.call(VCF_HEADER.format(self.path,
                                          'headers_' + self.path[-10:-7] + '.vcf'),
                        shell=True)
        subprocess.call(VCF_BODY.format(self.path,
                                        'body_' + self.path[-10:-7] + '.vcf'),
                        shell=True)
        subprocess.call(' '.join(['cat',
                                  'body_' + self.path[-10:-7] + '.vcf',
                                  '| java -jar cut.jar 1:9 >',
                                  'fixedf_' + self.path[-10:-7] + '.vcf']),
                        shell=True)

    def replace_calls(self):
        """
        Writes the new triplets of GL values for each SNP/variant
        in the vcf body.
        :return:
        """
        #self.samples_to_file()
        # transpose columns as lines
        subprocess.call(' '.join(['cat',
                                  'fixedf_' + self.path[-10:-7] + '.vcf',
                                  '| java -jar transpose.jar >',
                                  'T_fixedf_' + self.path[-10:-7] + '.vcf']),
                        shell=True)
        subprocess.call(' '.join(['cat',
                                  'new_data_' + self.kind + '_' + self.path[-10:-7] + '.vcf',
                                  '| java -jar transpose.jar >',
                                  'T_new_data_' + self.path[-10:-7] + '.vcf']),
                        shell=True)
        # append new data to the body
        subprocess.call(' '.join(['cat',
                                  'T_fixedf_' + self.path[-10:-7] + '.vcf',
                                  'T_new_data_' + self.path[-10:-7] + '.vcf',
                                  '>',
                                  'T_body_' + self.path[-10:-7] + '.vcf']),
                        shell=True)
        # re-transpose lines as columns
        subprocess.call(' '.join(['cat',
                                  'T_body_' + self.path[-10:-7] + '.vcf',
                                  '| java -jar transpose.jar >',
                                  'body_' + self.path[-10:-7] + '.vcf']),
                        shell=True)
        for pfx in ['new_data_' + self.kind + '_' ,
                    'T_new_data_',
                    'T_body_',
                    'T_fixedf_',
                    'fixedf_']:
            delete_file(pfx + self.path[-10:-7] + '.vcf')

    def reassemble_vcf(self):
        """
        Reassemble the VCF headers and and the new body.
        Delete intermediate vcf.
        :return:
        """
        subprocess.call(' '.join(['awk',
                                  '-v OFS="\t"',
                                  "'$1=$1'",
                                  'body_' + self.path[-10:-7] + '.vcf',
                                  '>',
                                  'body_tab_' + self.path[-10:-7] + '.vcf']),
                        shell=True)
        subprocess.call(' '.join(['cat',
                                  'headers_' + self.path[-10:-7] + '.vcf',
                                  'body_tab_' + self.path[-10:-7] + '.vcf',
                                  '>',
                                  'missing_' + self.kind + '_' + self.path[-10:-7] + '.vcf']),
                        shell=True) #'pooled_' + self.kind + '_' + self.path[-10:-7]
        # Gzip the output file
        #TODO: replace by bgzip/use bcftools to compress!
        _zip_vcf('missing_' + self.kind + '_' + self.path[-10:-7] + '.vcf')# 'pooled_' + self.kind + '_' + self.path[-10:-7]
        # # Index the output file for future processing
        # subprocess.call(' '.join(['bcftools',
        #                           'index -f',
        #                           '../data/pooled_GL_' + self.path[-10:-7] + '.vcf.gz']),
        #                 shell=True)

        print('{} bytes written'.format(file_size('pooled_'
                                                  + self.kind
                                                  + '_'
                                                  + self.path[-10:-7]
                                                  + '.vcf.gz')).ljust(80, '.'))

        for pfx in ['headers_', 'body_', 'body_tab_']:
            delete_file(pfx + self.path[-10:-7] + '.vcf')

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

parser = argparse.ArgumentParser(description='''Computes pooled GL/GT data for each subset.
 Subsets are written as new VCF files''')
parser.add_argument('start',
                    action='store',
                    help='ID of the first subset to be processed')
parser.add_argument('stop',
                    action='store',
                    help='ID of the last subset to be processed')
parser.add_argument('prefix',
                    action='store',
                    help='')
parser.add_argument('kind',
                    action='store',
                    help='gl/gt')
args = parser.parse_args()

files_list = [args.prefix + '{}.vcf.gz'.format(str(i).rjust(3, '0'))
              for i in range(int(args.start), int(args.stop) + 1)] # i in range(2, 159)
print('Processing the following files: ', files_list)

pll = mp.Pool(processes=4)
params = list(zip(files_list, [args.kind]*len(files_list)))
_ = pll.starmap(process_pool, params)

#TODO: specific processing for the remaining samples

