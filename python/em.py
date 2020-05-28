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
https://people.duke.edu/~ccc14/sta-663/EMAlgorithm.html
"""

xs = np.concatenate([np.tile((2, 2, 0), 4).reshape(4, 3),
                     np.tile((2, 1, 1), 12).reshape(12, 3),
                     np.tile((2, 0, 2), 4).reshape(4, 3),
                     np.tile((1, 3, 0), 4).reshape(4, 3),
                     np.tile((1, 2, 1), 12).reshape(12, 3),
                     np.tile((1, 1, 2), 12).reshape(12, 3),
                     np.tile((1, 0, 3), 4).reshape(4, 3),
                     np.tile((0, 4, 0), 1).reshape(1, 3),
                     np.tile((0, 0, 4), 1).reshape(1, 3),
                     np.tile((0, 2, 2), 6).reshape(6, 3),
                     np.tile((0, 1, 3), 4).reshape(4, 3),
                     np.tile((0, 3, 1), 4).reshape(4, 3)
                     ])

valid = np.concatenate([np.repeat(2/4, 4),
                        np.repeat(4/12, 12),
                        np.repeat(2/4, 4),
                        np.repeat(1, 4),
                        np.repeat(1, 12),
                        np.repeat(1, 12),
                        np.repeat(1, 4),
                        np.repeat(1, 1),
                        np.repeat(1, 1),
                        np.repeat(1, 6),
                        np.repeat(1, 4),
                        np.repeat(1, 4)
                        ])

# theta trinomial distribution
thetas = {'hwe': np.array([0.25, 0.5, 0.25]),
          'uni': np.array([0.333, 0.333, 0.333])
          }

weights = {'hwe': np.array([1.0, 2.0, 1.0]),
           'uni': np.array([1.0, 1.0, 1.0])
           }

# data observed for the unknown genotypes over all possible completions
counts = np.sum(xs, axis=0)
# hence the likelihood per sample (layout)
lkh = np.apply_along_axis(lambda x: x/np.sum(x, axis=0), axis=1, arr=xs)
# proportions observed, the proba distribution behind is unknown
probs = counts / np.sum(counts)

# expected counts if the distribution was theta: expectation formula
Pr_cond = {}
E_counts = {}
proportions = {}
posteriors = {'hwe': [],
              'uni': []
              }
priors = {'hwe': [thetas['hwe']],
          'uni': [thetas['uni']]
          }

max_iter = 30

pdic = {'hwe': thetas['hwe'], 'uni': thetas['uni']}
for i in range(max_iter):
    for k, v in pdic.items():
        # E step
        #Pr_cond[k] proba cond of each layout given its frequency validity and the prior
        Pr_cond[k] = np.multiply(np.multiply.reduce(np.power(v, xs), axis=1),
                                 valid)
                                 # np.repeat(valid, 3).reshape((len(valid), 3)))
        # proportion of each layout in the mix resulting in the pattern
        proportions[k] = np.apply_along_axis(lambda x: x / np.sum(x), axis=0, arr=Pr_cond[k])

        # M step: resize each layout by proportion and count genotypes over all layouts
        E_counts[k] = np.multiply(np.repeat(proportions[k], 3).reshape((len(xs), 3)),
                                  xs).sum(axis=0)
        E_freq = E_counts[k] / E_counts[k].sum()

        E_w = np.multiply(weights[k],
                          np.divide(E_freq, v))  # why?
        pr_w = np.multiply(1/weights[k], E_counts[k])

        freq_up = np.apply_along_axis(lambda x: x/np.sum(x), axis=0, arr=E_w)
        freq_down = np.apply_along_axis(lambda x: x/np.sum(x), axis=0, arr=pr_w)

        new_prior = E_w / np.sum(E_w)  #(1 / len(xs)) * np.sum(freq_up, axis=0)
        posterior = pr_w / np.sum(pr_w)  #(1 / len(xs)) * np.sum(freq_down, axis=0)
        priors[k].append(new_prior)
        pdic[k] = priors[k][-1]
        posteriors[k].append(posterior)
