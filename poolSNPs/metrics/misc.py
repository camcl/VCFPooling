"""
Shannon's diversity index for assesssing genotypes spreading/uncertainty about the most likely genotype
(overall within a population or locally for one individual.
For instance, the case GL=[0.33, 0.33, 0.33] stands for the situation where the diversity is maximum.

Implement metrics for measuring diversity in population genotypes: Mean Squre diff? at population level
"""

import numpy as np
from typing import *
from collections import Counter
import os, sys

rootdir = os.path.dirname(os.path.dirname(os.getcwd()))
sys.path.insert(0, rootdir)

from VCFPooling.poolSNPs import dataframe as vcfdf
from VCFPooling.persotools.files import *

ArrayLike = NewType('ArrayLike', Union[Sequence, List, Set, Tuple, Iterable, np.ndarray, int, float, str])


def normalize(v: np.ndarray) -> np.ndarray:
    normv = np.linalg.norm(v)
    if normv == 0.0:
        return v
    else:
        return v/normv

def shannons_index(a: np.ndarray) -> float:
    """
    see scipy.stats.entropy
    """
    hfunc = np.vectorize(lambda x: -x * np.log2(x) if x != 0.0 else 0)
    h = np.sum(hfunc(a), axis=0)
    return h


def gt_diversity(a: np.ndarray) -> float:
    """
    Diversity in GT at a given marker. GT trinary encoded and possibly missing.
    """
    scores = Counter(a)
    # print(scores.items())
    p = [s/len(a) for s in scores.values()]
    return shannons_index(p)


def gl_diversity(a: np.ndarray) -> float:
    """
    Diversity in GL at a given marker.
    a has shape (n, 3) where is the number of samples
    """
    p = np.mean(a, axis=-1)  # 'mean' GL for the pop
    return shannons_index(p)


class Diversity(object):
    """

    """
    def __init__(self, filepath: FilePath, format: str = None, idx: str = 'id'):
        self.obj = vcfdf.PandasMixedVCF(filepath, format=format, indextype=idx)
        self.fmt = format

    def markers_diversity(self):
        if self.fmt == 'GT':
            df = self.obj.trinary_encoding()
            dvs = df.apply(gt_diversity, axis=1)
            return dvs
        if self.fmt == 'GL':
            pass


if __name__ == '__main__':
    bas10000 = 'IMP.chr20.snps.gt.chunk1000.vcf.gz'
    path = '/home/camille/1000Genomes/data/gt/stratified/'
    div0 = Diversity(path + bas10000, format='GT').markers_diversity()

    gt10000 = 'IMP.chr20.pooled.snps.gt.chunk1000.vcf.gz'
    gtpath = '/home/camille/1000Genomes/data/gt/'
    div1 = Diversity(gtpath + gt10000, format='GT').markers_diversity()

    imp10000 = 'IMP.chr20.pooled.imputed.gt.chunk10000.vcf.gz'
    glpath = '/home/camille/1000Genomes/data/gl/gl_adaptive/all_snps_all_samples/'
    div2 = Diversity(glpath + imp10000, format='GT').markers_diversity()
