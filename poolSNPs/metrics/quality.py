"""
Metrics for assessing imputation quality?

Het/Hom ratio

Improving imputation quality in BEAGLE for crop and
livestock data
    - switch-error rate for imputation quality
    - idem ?

Comparison and assessment of family-
and population-based genotype imputation
methods in large pedigrees
    - Mean squared correlation (R^2) = Pearson’s squared correlation: [0, 1].
    - concordance rate (CR): overestimates the imputation accuracy for rare variants
    - imputation quality score (IQS): agreement ratio, (-Inf, 1],
    based on the Kappa statistic


A New Statistic to Evaluate Imputation Reliability
IQS:
The computation of IQS requires the posterior probabilities of AA, AB and BB as output by the imputation program.
--> with Beagle 4.1: gprobs, in dic['corr'] files as GT:DS:GP ex. 0|0:0.51:0.55,0.38,0.07
gprobs=[true/false]specifies  whether  a  GP  (genotype  probability)
format  field  will  be included in the output VCF file (default: gprobs=true)


MaCH: Using Sequence and Genotype Data to Estimate Haplotypes and Unobserved Genotypes
1 SNP with alleles A and B. Let n_A/A , n_A/B , n_B/B = number of times
each possible genotype was sampled after I = n_A/A + n_A/B + n_B/B iterations
Most likely genotype = genotype that was sampled most frequently
Expected number of counts of allele A: g = (2*n_A/A + n_A/B)/I

    1) Genotype Quality Score: GQS = n_IG /I,
    n_IG = number of iterations where the given genotype was selected as the most likely one
    This quantity can be averaged over all genotypes for a
    particular marker to quantify the average accuracy of imputation for that marker

    2) Accuracy: alpha = sum(GQS_i, i=1:N)/N,
    N number of individuals

    3) R²: E(r² with true genotypes) = Var(g)/((4*n_A/A + n_A/B)/I - [(2*n_A/A + n_A/B)/I]²),
    Estimated r² with true genotypes, Var(g) be the variance of estimated genotype:
    a better measure of imputation quality for a marker is the estimated r² between
    true allele counts and estimated allele counts. This quantity can be estimated by
    comparing the variance of the estimated genotype scores with what would be expected if
    genotype scores were observed without error.


https://en.wikipedia.org/wiki/Cohen's_kappa
Cohen's kappa coefficient K is a statistic which measures inter-rater agreement for qualitative (categorical) items.
Generally thought to be more robust than simple percent agreement calculation, as that coeff takes into account
the possibility of agreement occuring by chance. But: difficult to interpret indices of agreement?
If no agreement between the raters other than the one that would be expected by chance: K = 0,
If K < 0: there is no effective agreement between the raters or the agreement is worse than random.
Weighted kappa K exists.
κ's tendency to take the observed categories' frequencies as givens, which can make it unreliable for measuring
agreement in situations such as the diagnosis of rare diseases. In these situations, κ tends to underestimate
the agreement on the rare category.


https://scikit-learn.org/stable/modules/generated/sklearn.metrics.cohen_kappa_score.html#sklearn.metrics.cohen_kappa_score
Implementation of Cohen's kappa:
sklearn.metrics.cohen_kappa_score(y1, y2, labels=None, weights=None, sample_weight=None)
y1, y2: arrays of same length (n_samples)


http://courses.washington.edu/cmling/lab7.html
Using the python interpreter and the nltk metrics package, calculate inter-annotator agreement (both kappa and alpha).
Note that AnnotationTask is a type of object, with methods kappa() and alpha().
When you call nltk.metrics.AnnotationTask() it returns an object of that type, which in the example below is stored
in the variable task. See: http://www.nltk.org/api/nltk.metrics.html

import nltk
toy_data = [
['1', 5723, 'ORG'],
['2', 5723, 'ORG'],
['1', 55829, 'LOC'],
['2', 55829, 'LOC'],
['1', 259742, 'PER'],
['2', 259742, 'LOC'],
['1', 269340, 'PER'],
['2', 269340, 'LOC']
]
task = nltk.metrics.agreement.AnnotationTask(data=toy_data)
task.kappa()
task.alpha()

The nltk metrics package also provides for calculating and printing confusion matrices, a way of displaying which labels
 were 'mistaken' for which other ones. Unfortunately, this functionality requires a different format for the input.
 In particular, it wants two lists of labels (in the same order).
"""

import os
import numpy as np
from scipy.stats import pearsonr
from typing import *

from scripts.poolSNPs import parameters as prm
from scripts.poolSNPs import patterns as pat
from scripts.poolSNPs.alleles import alleles_tools as alltls
from persotools.debugging import *
from persotools.files import *

ArrayLike = NewType('ArrayLike', Union[Sequence, List, Set, Tuple, Iterable, np.ndarray, int, float, str])


def correlation_r2(vcf_true: FilePath, vcf_imputed: FilePath) -> Tuple:
    try:
        x = alltls.get_aaf(vcf_true, id='chrom:pos')
        y = alltls.get_aaf(vcf_imputed, id='chrom:pos')
    except IOError:
        print('VCF files must provide GT-formatted genotypes')
        x, y = None, None

    r, p_value = pearsonr(x, y)
    return r, p_value

#TODO: implement method for extracting GP field