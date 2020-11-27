import numpy as np
import itertools
from scipy.linalg import block_diag

import matplotlib.pyplot as plt

"""
Playground for self-consistent genotype decoding population-wide
"""


# Design matrix

pools_size = 4
pools_nb = 8
block_size = 16
blocks_nb = 1
ndim = 2

m = np.zeros((pools_nb, block_size), dtype=int)
for i in range(int(pools_nb/ndim)):
    j = i * pools_size
    m[i, j:j+pools_size] = [1]*pools_size
for i in range(int(pools_nb/ndim), pools_nb):
    j = i - pools_size
    m[i, [j+k*pools_size for k in range(pools_size)]] = 1
b = itertools.repeat(m, blocks_nb)
M = block_diag(*b)


# Cross-entropy calculation

prior = np.array([0.25, 0.5, 0.25])

mylog12 = np.vectorize(lambda x: np.log(x) if x > 1e-12 else np.log(1e-12))

prodlog = np.vectorize(lambda x, y: -np.multiply(np.asarray(x, dtype=float),
                                                 np.log(np.asarray(y, dtype=float))
                                                 )
                       )

myprodlog = np.vectorize(lambda x, y: -np.multiply(np.asarray(x, dtype=float),
                                                   mylog12(np.asarray(y, dtype=float))
                                                   )
                         )


# Genotypes (GP format)

nb_samples = block_size * blocks_nb

hom_ref_only = np.tile([1.0, 0.0, 0.0], nb_samples).reshape(nb_samples, 3)
hom_ref_but_one = np.vstack([hom_ref_only[:-1, :], np.array([0.0, 1.0, 0.0])])

rd = np.random.random((nb_samples, 3))
mixed = rd / rd.sum(axis=-1)[:, np.newaxis]

a = hom_ref_but_one

#  Encoding

scores_p = np.tensordot(M, a, axes=1)
genos_p = scores_p / scores_p.sum(axis=-1)[:, np.newaxis]


#  Decoding

scores_i = np.tensordot(M.T, genos_p, axes=1)
genos_i = scores_i / scores_i.sum(axis=-1)[:, np.newaxis]

posterior = genos_i.sum(0) / genos_i.sum()


# Cross-entropy?
E0 = myprodlog(a[:, 0], genos_i[:, 0]).sum()
E1 = myprodlog(a[:, 1], genos_i[:, 1]).sum()
E2 = myprodlog(a[:, 2], genos_i[:, 2]).sum()

entropy = np.array([E0, E1, E2])
# /1\ log(0) undefined issues

posterior = np.power(prior, entropy)

# Signal denoising?

import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import OrthogonalMatchingPursuit
from sklearn.linear_model import OrthogonalMatchingPursuitCV

n_nonzero_coefs = [np.round(a[:, 0].sum()).astype(int),
                   np.round(a[:, 1].sum()).astype(int),
                   np.round(a[:, 2].sum()).astype(int)
                   ]

# omp = OrthogonalMatchingPursuit(n_nonzero_coefs=n_nonzero_coefs[0])
# omp.fit(a[:, 0].reshape(-1, 1), genos_i[:, 0].reshape(-1, 1))
# coef = omp.coef_
# idx_r, = coef.nonzero()


# Image representation

plt.imshow(a)
plt.show()

plt.imshow(genos_i)
plt.show()

from skimage.restoration import (denoise_wavelet, estimate_sigma)
from skimage import data, img_as_float
from skimage.util import random_noise
from skimage.metrics import peak_signal_noise_ratio
from skimage.color import rgb2hsv

import cv2

hsv_genos = rgb2hsv(genos_i)
plt.imshow(hsv_genos)
plt.show()

s = hsv_genos[:, 1]

img = genos_i
percent = 0.5
#percent = 0.25
#percent = 0

# desaturate
s_desat = cv2.multiply(s, percent).astype(np.uint8)
