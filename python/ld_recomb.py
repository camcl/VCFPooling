import itertools as it
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
plt.style.use('ggplot')
#np.random.seed(1234)

# np.set_printoptions(formatter={'all':lambda x: '%.3f' % x})

from numpy.core.umath_tests import matrix_multiply as mm
from scipy.optimize import minimize
from scipy.stats import bernoulli, binom

"""
http://stephenslab.uchicago.edu/assets/papers/Li2003.pdf
"""

rho = np.arange(0.01, 0.72, 0.02)  # recombination rate
d = np.arange(0.01, 100.1, 0.1)  # physical distance beween markers
k = np.arange(10, 101, 10)  # number of template haplotypes


def pr_transition(x1, x2, recomb_rate=0.3, dist=10, nb_haplo=10):
    if x1 == x2:
        pr_trans = np.exp(-recomb_rate * dist / nb_haplo) + (1/nb_haplo) * (1 - np.exp(-recomb_rate * dist / nb_haplo))
    else:
        pr_trans = (1/nb_haplo) * (1 - np.exp(-recomb_rate * dist / nb_haplo))

    return pr_trans


rho_var_stay = np.vectorize(lambda x: pr_transition(1, 1, recomb_rate=x))
rho_var_move = np.vectorize(lambda x: pr_transition(1, 2, recomb_rate=x))

d_var_stay = np.vectorize(lambda x: pr_transition(1, 1, dist=x))
d_var_move = np.vectorize(lambda x: pr_transition(1, 2, dist=x))

k_var_stay = np.vectorize(lambda x: pr_transition(1, 1, nb_haplo=x))
k_var_move = np.vectorize(lambda x: pr_transition(1, 2, nb_haplo=x))

fig, ax = plt.subplots(1, 3, sharey=True)
ax[0].scatter(rho, rho_var_stay(rho), c='b')
ax[0].scatter(rho, rho_var_move(rho), c='g')
ax[0].set_ylim(0.0, 1.1)
ax[1].scatter(d, d_var_stay(d), c='b')
ax[1].scatter(d, d_var_move(d), c='g')
ax[1].set_ylim(0.0, 1.1)
ax[2].scatter(k, k_var_stay(k), c='b')
ax[2].scatter(k, k_var_move(k), c='g')
ax[2].set_ylim(0.0, 1.1)
plt.show()
