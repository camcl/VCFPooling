import sys, os
import argparse
import timeit

# force PYTHONPATH to look into the project directory for modules
bin_dir = os.path.dirname(os.getcwd())
print(bin_dir)
sys.path.insert(0, bin_dir)

import poolSNPs.parameters as prm
from poolSNPs import poolvcf
from poolSNPs import pybcf

'''
Applies pooling simulation to a VCF file

* the input VCF file contains variants of type SNP only, with genotype formatted as GT,
* the output VCF file contains the same variants, formatted as GL,
* the number of samples is a multiple of 16 (block's size in the DNA Sudoku design implemented).
The decoding step of pooling is adaptive with respect to the pooling pattern observed (see README.md)
* the samples are assumed to be sorted in row-order flattened blocks order e.g. the 16 first columns in the VCF file
correspond  to the samples assigned to the first block. 
Samples 1-4 form the first pool in the block, samples 5-8 the second pool, and so on.
* the output file has to be written to unbgzipped format (.vcf) and then compressed to 
bgzipped format (.vcf.gz) with bcftools.

For VCF-file bigger than some dozen of thousands of variants, pooling can be parallelized.

Command line usage (assuming the current directory is VCFPooling/examples
$ python3 -u pooling-ex.py <path-to-file-in> <path-to-file-out>
'''

### COMMAND-LINE PARSING AND PARAMETERS
parser = argparse.ArgumentParser(description='Run pooling simulation'
                                             'on the whole set of samples')
parser.add_argument('pathin', metavar='in', type=str, help='File to pool', default=None)
parser.add_argument('pathout', metavar='out', type=str, help='Pooled file', default=None)

argsin = parser.parse_args()
filin = argsin.pathin
filout = argsin.pathout
plookup = os.path.join(os.getcwd(), 'adaptive_gls.csv')  # look-up table for converting pooled GT to GL

print('\n'.ljust(80, '*'))
print('Input file = {}'.format(os.path.expanduser(argsin.pathin)))
print('Output file = {}'.format(os.path.expanduser(argsin.pathout)))
print('\n'.rjust(80, '*'))

# make sure to write to .vcf
if filout.endswith('.gz'):
    vcfout = filout[:-3]

### SIMULATE POOLING
start = timeit.default_timer()
poolvcf.pysam_pooler(filin, vcfout, plookup, os.getcwd())

print('\r\nTime elapsed --> ', timeit.default_timer() - start)
