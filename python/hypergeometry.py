import numpy as np
import scipy.stats as stats
import pandas as pd
import matplotlib.pyplot as plt
from itertools import starmap
from functools import partial
from cyvcf2 import VCF

from scripts.VCFPooling.poolSNPs.alleles import alleles_tools as alltls
from scripts.VCFPooling.poolSNPs import parameters as prm
from persotools.files import *

"""
Plot AAF after pooling (pooled)
and AAF in the original dataset (theoretical).
Compare with hypergeometrical cumulative density function

Ref article: Identifying Rare Variants With Optimal Depth of Coverage and Cost-Effective Overlapping Pool Sequencing.
Chang-Chang Cao, Cheng Li, Zheng Huang, Xin Ma, and Xiao Sun. 2013.

For a population of n samples dividing into B blocks with
each has n B samples, suppose the frequency of a rare variant is
p v , the number of positive samples (d B ) in each block follows
a binomial distribution B(n B , p v ) and the probability p B that
the positive samples are no more than d B is:
When p v is unknown but there exist at most d positive
samples in the population, the number of positive samples in
each block follows a hypergeometric distribution H(n B , d, n)
and the probability p B that the positive samples are no more
than d B is
"""

N = 2496  # number of samples
B = 156  # nb of blocks
nB = 16  # nb of idv per block
aaf = np.arange(0.0, 1.05, 0.05)
maf = np.arange(0.0, 0.55, 0.05)
lwf = np.arange(0.0, 0.1, 0.005)
carrier = np.vectorize(lambda x: 2 * x * (1-x) + x**2)
acf = carrier(aaf)  # alternate allele carrier frequency
mcf = carrier(maf)  # alternate allele carrier frequency
lcf = carrier(lwf)  # carriersat low frequencies
k = 1  # decoding power, number of carriers per pool

# proba there are 2 or more carriers a pool given AAF
# hypergeom.sf(k, N, N*acf, nB, loc=0)
print('Tot number of carriers: ', np.round(N*acf))
sf = stats.hypergeom.sf(k, N, np.round(N*acf), nB, loc=nB//2)
sfsym = stats.hypergeom.sf(nB-k+1, N, np.round(N*acf), nB, loc=nB//2)  # symmetry of ref/alt
vecsf = np.subtract(sfsym, sf)
print(vecsf)
cdf = stats.hypergeom.cdf(k, N, np.round(N*acf), nB, loc=nB//2)
cdfsym = stats.hypergeom.cdf(nB-k+1, N, np.round(N*acf), nB, loc=nB//2)  # symmetry of ref/alt
veccdf = np.subtract(cdfsym, cdf)

cumul = np.add(cdf, sfsym)
probs = cumul
tweak = np.multiply(aaf, probs)
print('tweaked AAF', tweak)

plt.plot(acf, cumul)
plt.xlabel("Alternate allele carrier frequency")
plt.ylabel("Cumulated probabilty encountering x% alternate allele carriers")
plt.title("Cumulative Distribution Function of a hypergeometric probability law")
plt.show()

fig, ax = plt.subplots(2, 1)
ax[0].plot(aaf, cumul)
ax[1].plot(aaf, tweak)
plt.show()
