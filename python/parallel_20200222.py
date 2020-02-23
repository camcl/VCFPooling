import sys, os
from typing import *
import numpy as np
from cyvcf2 import VCF
from itertools import starmap, repeat
import shutil
import multiprocessing as mp
import argparse

home_dir = os.path.expanduser("~")
proj_dir = os.path.join(home_dir, '1000Genomes')
sys.path.insert(0, proj_dir)

from scripts.VCFPooling.poolSNPs import parameters as prm
from scripts.VCFPooling.poolSNPs import beagle_tools as bgltls
from scripts.VCFPooling.poolSNPs import pybcf
from scripts.VCFPooling.poolSNPs import chunkvcf as chkf
from scripts.VCFPooling.poolSNPs.alleles import alleles_tools as alltls

from persotools.files import delete_file, mkdir
from persotools.struct import NamedDict

'''
Parallelized file processing: Read main vcf and write chunks

Steps:
* 
*
*
'''




if __name__ == '__main__':
    os.chdir('/home/camille/1000Genomes/data/gl/gl_adaptive/all_snps_all_samples')

    gtframe = alltls.PandasMixedVCF('IMP.chr20.pooled.beagle2.gl.chunk10000.corr.vcf.gz', format='GP')
    varonly = chkf.VariantCallGenerator('IMP.chr20.pooled.beagle2.gl.chunk10000.corr.vcf.gz', format='GP')
    varchunk = chkf.VariantChunkGenerator('IMP.chr20.pooled.beagle2.gl.chunk10000.corr.vcf.gz',
                                            format='GP',
                                            chunksize=1000)

    chunkpack = varchunk.chunkpacker()
    for i, chk in enumerate(chunkpack):
        print(i, len([*chk]))
