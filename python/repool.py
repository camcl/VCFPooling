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
blocks_nb = 10
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

nb_patterns = 3

hom_ref_only = np.tile([1.0, 0.0, 0.0], block_size).reshape(block_size, 3)  # pattern 1
hom_ref_but_one = np.vstack([hom_ref_only[:-1, :], np.array([0.0, 1.0, 0.0])])  # pattern 2
hom_ref_2_alt = np.vstack([np.array([0.0, 0.0, 1.0]), hom_ref_only[1:-1, :], np.array([0.0, 1.0, 0.0])])  # pattern 3

rd = np.random.random((nb_samples, 3))
mixed = rd / rd.sum(axis=-1)[:, np.newaxis]

a = np.zeros((blocks_nb, block_size, 3))
block_choice = np.random.binomial(nb_patterns - 1, 0.4, size=blocks_nb)

a[block_choice.flatten() == 1] = hom_ref_but_one
a[block_choice.flatten() == 0] = hom_ref_only
a[block_choice.flatten() == 2] = hom_ref_2_alt

a = a.reshape(nb_samples, 3)

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
# /!\ log(0) undefined issues

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

plt.rcParams["image.cmap"] = 'binary'
plt.rcParams["figure.figsize"] = (2 * blocks_nb, 2)

plt.imshow(a.T)
plt.show()

plt.imshow(genos_i.T)
plt.show()

import pybeads as be
from skimage.color import rgb2hsv

import cv2

# hsv_genos = rgb2hsv(genos_i)
# plt.imshow(hsv_genos)
# plt.show()

y = 1 - genos_i[:, 0]


def sigmoid(x):
    return 1 / (1 + np.exp(-x))


# y is your data

xscale_l, xscale_r = 100, 100
dx = 1
y_l = y[0]*sigmoid(1/xscale_l*np.arange(-5*xscale_l, 5*xscale_l, dx))
y_r = y[-1]*sigmoid(-1/xscale_r*np.arange(-5*xscale_r, 5*xscale_r, dx))
y_ext = np.hstack([y_l, y, y_r])

# try changing fc and lam0-2, amp if these dont' fit your data
fc = 0.006
d = 1
r = 6
amp = 0.8
lam0 = 0.5 * amp
lam1 = 5 * amp
lam2 = 4 * amp
Nit = 15
pen = 'L1_v2'

# signal_est, bg_est, cost = be.beads(y_ext, d, fc, r, Nit, lam0, lam1, lam2, pen, conv=None)
