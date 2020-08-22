"""
Metrics for assessing imputation quality

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
    This ity can be averaged over all genotypes for a
    particular marker to ify the average accuracy of imputation for that marker

    2) Accuracy: alpha = sum(GQS_i, i=1:N)/N,
    N number of individuals

    3) R²: E(r² with true genotypes) = Var(g)/((4*n_A/A + n_A/B)/I - [(2*n_A/A + n_A/B)/I]²),
    Estimated r² with true genotypes, Var(g) be the variance of estimated genotype:
    a better measure of imputation quality for a marker is the estimated r² between
    true allele counts and estimated allele counts. This ity can be estimated by
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

import os, sys
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, zscore
from scipy.special import softmax
from sklearn import metrics, preprocessing
from typing import *

rootdir = os.path.dirname(os.path.dirname(os.path.dirname(os.getcwd())))
sys.path.insert(0, rootdir)

from VCFPooling.poolSNPs import dataframe as vcfdf
from VCFPooling.poolSNPs.metrics.misc import normalize, min_max_scale
from VCFPooling.persotools.files import *

ArrayLike = NewType('ArrayLike', Union[Sequence, List, Set, Tuple, Iterable, np.ndarray, int, float, str])

#TODO: evaluate phase/switch rate


class QualityGT(object):
    """
    Implement different methods for assessing imputation performance:
    * accuracy and recall per variant per genotype (cross-table)
    * correlation per variant and/or per sample between imputed and true genotypes
    * difference per variant and/or per sample between imputed and true genotypes
    * allele dosage
    """
    def __init__(self, truefile: FilePath, imputedfile: FilePath, ax: object, idx: str = 'id'):
        self.trueobj = vcfdf.PandasMixedVCF(truefile, format='GT', indextype=idx)
        self.imputedobj = vcfdf.PandasMixedVCF(imputedfile, format='GT', indextype=idx)
        self._axis = ax
        #TODO: index properties and verification

    @property
    def axis(self):
        return self._axis

    @axis.setter
    def set_axis(self, ax):
        if ax == 0 or ax == 'variants':
            self._axis = 0
        elif ax == 1 or ax == 'samples':
            self._axis = 1
        else:
            self._axis = None

    @staticmethod
    def square(x):
        return x ** 2

    def pearsoncorrelation(self) -> pd.Series:
        """
        Compute Pearson's correlation coefficient between true and imputed genotypes.
        Correlation between variants (ax=1 i.e. mean genotypes along samples axis),
        or correlation between samples (ax=0 i.e. mean genotypes along variant axis),
        or global correlation (ax=None i.e. mean of flattened array)
        :return: correlation coefficients and p-value for each
        """
        #TODO: replace by Allele Frequency correlation as described in Beagle09?
        true = self.trueobj.trinary_encoding().values
        imputed = self.imputedobj.trinary_encoding().values
        scorer = lambda t: pearsonr(t[0], t[1])[0]  # keeps only correlation, not p-value
        score = list(map(scorer, zip(true, imputed)))
        # astype(str) casts Series type to discrete classes
        rsqr = pd.Series(score, index=self.trueobj.variants, name='r_squared').apply(self.square)
        # squared correlation
        return rsqr

    def diff(self) -> pd.DataFrame:
        """
        Compute absolute genotype difference element-wise i.i per variant per sample
        :return: absolute difference true vs. imputed genotypes
        """
        truedf = self.trueobj.trinary_encoding()
        imputeddf = self.imputedobj.trinary_encoding()
        absdiffdf = truedf.sub(imputeddf).abs()
        return absdiffdf

    def concordance(self) -> pd.Series:
        """
        Compute concordance between true and imputed genotypes
        i.e. 1 - the Z-norm of the absolute difference of true vs. imputed genotypes?
        :return:
        """
        # absdiff = self.diff()  # equals 0 when true = imputed, else can be 1 or 2 (very unlikely 2?)
        absdiff = self.diff() / 2  # restricts values to 0.0, 0.5, 1.0
        # absdiffnorm = absdiff.apply(min_max_scale, axis=1, raw=True)  # homemade minmax scaler
        absdiffnorm = absdiff.apply(preprocessing.minmax_scale, axis=1, raw=True)  # sklearn minmax scaler
        discord_score = absdiffnorm.mean(axis=1)  # discordance
        concord_score = 1 - discord_score  # concordance = 1 - discordance
        concord = pd.Series(concord_score, index=self.trueobj.variants, name='concordance')
        return concord

    @staticmethod
    def expectation(a: np.ndarray, freq: np.ndarray):
        """
        sum(Pr(G=x)*x)
        :param a:
        :param freq:
        :return:
        """
        return np.multiply(a, freq).sum()

    def alleledosage(self) -> Tuple[pd.Series]:
        # TODO: add  by Standardized Allele Frequency Error as described in Beagle09?
        """
        Compute mean genotype along the given axis. Equivalent to AAF computation.
        Makes sense only accross a population i.e. mean values along samples axis.
        Allele dosage = 2 * AAF, for a diploid organism
        :return:
        """
        truedos = self.trueobj.trinary_encoding().values.mean(axis=1)
        imputeddos = self.imputedobj.trinary_encoding().values.mean(axis=1)
        strue = pd.Series(truedos, index=self.trueobj.variants, name='truedos')
        simputed = pd.Series(imputeddos, index=self.imputedobj.variants, name='imputeddos')

        return strue, simputed

    @property
    def precision(self, avg: str = 'weighted') -> pd.Series:
        """
        Compute precision score for the imputed genotypes.
        The precision is the ratio tp / (tp + fp) where tp is the number of true positives and
        fp the number of false positives. The precision is intuitively the ability of the classifier
        not to label as positive a sample that is negative. The best value is 1 and the worst value is 0.
        :param avg: 'weighted' needed for multiclass classification
        :return:
        """
        true = self.trueobj.trinary_encoding().values
        imputed = self.imputedobj.trinary_encoding().values
        scorer = lambda t: metrics.precision_score(t[0].astype(str),
                                                   t[1].astype(str),
                                                   average=avg)
        score = list(map(scorer, zip(true, imputed)))
        # astype(str) casts Series type to discrete classes
        return pd.Series(score, index=self.trueobj.variants, name='precision_score')

    @property
    def accuracy(self) -> pd.Series:
        """
        Compute accuracy score for the imputed genotypes.
        In multilabel classification, this function computes subset accuracy i.e. the number of exact true matches.
        The accuracy is the ratio tp / (tp + fp + tn + fn) for each class.
        Equal to Jaccard index in the case of multilabel classification tasks.
        Jaccard similarity coefficient is defined as the size of the intersection
        divided by the size of the union of two label sets.
        :return:
        """
        true = self.trueobj.trinary_encoding().values
        imputed = self.imputedobj.trinary_encoding().values
        scorer = lambda t: metrics.accuracy_score(t[0].astype(str),
                                                  t[1].astype(str))
        score = list(map(scorer, zip(true, imputed)))
        # astype(str) casts Series type to discrete classes
        return pd.Series(score, index=self.trueobj.variants, name='accuracy_score')

    @property
    def recall(self, avg: str = 'weighted') -> pd.Series:
        """
        Compute recall score for the imputed genotypes.
        The recall is the ratio tp / (tp + fn) where tp is the number of true positives and
        fn the number of false negatives. The recall is intuitively the ability of the classifier
        to find all the positive samples.
        The best value is 1 and the worst value is 0.
        :return:
        """
        true = self.trueobj.trinary_encoding().values
        imputed = self.imputedobj.trinary_encoding().values
        scorer = lambda t: metrics.recall_score(t[0].astype(str),
                                                t[1].astype(str),
                                                average=avg)
        score = list(map(scorer, zip(true, imputed)))
        # astype(str) casts Series type to discrete classes
        return pd.Series(score, index=self.trueobj.variants, name='recall_score')

    @property
    def f1_score(self, avg: str = 'weighted') -> pd.Series:
        """
        F1-score for the genotypes
        :return:
        """
        true = self.trueobj.trinary_encoding().values
        imputed = self.imputedobj.trinary_encoding().values
        scorer = lambda t: metrics.f1_score(t[0].astype(str),
                                            t[1].astype(str),
                                            average=avg)
        score = list(map(scorer, zip(true, imputed)))
        # astype(str) casts Series type to discrete classes
        return pd.Series(score, index=self.trueobj.variants, name='f1_score')


class QualityGL(object):
    """
    Implement cross-entropy method for assessing imputation performance from GL
    """
    def __init__(self, truefile: FilePath, imputedfile: FilePath, ax: object, fmt: str = 'GP', idx: str = 'id'):
        self.trueobj = vcfdf.PandasMixedVCF(truefile, format='GL', indextype=idx)
        self.imputedobj = vcfdf.PandasMixedVCF(imputedfile, format=fmt, indextype=idx)
        self._axis = ax
        #TODO: index properties and verification

    @property
    def axis(self):
        return self._axis

    @axis.setter
    def set_axis(self, ax):
        if ax == 0 or ax == 'variants':
            self._axis = 0
        elif ax == 1 or ax == 'samples':
            self._axis = 1
        else:
            self._axis = None

    def logfill(self, x):
        """
        Adjust values for cross-entropy calculation
        Ex. rs1836444  -0.0  inf  NaN -> rs1836444  0.0  5.0  0.0
        """
        if x == -0.0 or x == np.nan:
            return 0.0
        elif x == np.inf:
            return 5.0
        else:
            return x

    def intergl_entropy(self, g_true: pd.Series, g_pred: pd.Series) -> pd.Series:
        """
        Compute entropy from two GL series for a sample as
        E = -sum(p_true * log(p_imputed), sum over the 3 GL values at every mmarker
        p_imputed set to 10^-5 if equal to 0.0
        """
        g_true = pd.Series(g_true)
        g_pred = pd.Series(g_pred)  # comes as tuples of str
        dftrue = pd.DataFrame.from_records(g_true.values,
                                           index=g_true.index,
                                           columns=['RR', 'RA', 'AA']).astype(float)
        dfpred = pd.DataFrame.from_records(g_pred.values,
                                           index=g_pred.index,
                                           columns=['RR', 'RA', 'AA']).astype(float)
        prodlog = lambda x, y: -np.multiply(np.asarray(x, dtype=float),
                                            np.log(np.asarray(y, dtype=float))
                                            )
        # GP are not logged
        g_entro = dftrue.combine(dfpred, prodlog).applymap(self.logfill).fillna(0.0)
        return g_entro.sum(axis=1).rename(g_true.name)

    @property
    def cross_entropy(self) -> pd.Series:
        """
        For genotypes likelihoods
        Entropy for the genotypes, aCROSS two populations.
        Not confuse with intrapop entropy
        entropy = alpha * sum(p_true * log(p_imputed) for every GL for every sample) at 1 marker
        :return:
        """
        true = self.trueobj.genotypes()
        imputed = self.imputedobj.genotypes()
        entro = true.combine(imputed, self.intergl_entropy)
        score = entro.mean(axis=1)
        return pd.Series(score, index=self.trueobj.variants, name='cross_entropy')


class QuantilesDataFrame(object):
    """
    Builds a DataFrame with quantiles values over bins for a given Quality data Series object
    and its AF_INFO/MAF_INFO data Series
    """
    def __init__(self, dX: pd.DataFrame, dY: pd.Series, bins_step: float = 0.01):
        # assert dX.columns[0] in ['maf_info', 'af_info', 'maf', 'aaf']
        self.dX = dX
        self.dY = dY.to_frame()
        self.x_data = dX.columns[0]
        self.y_data = dY.name

        self.y_bins = np.arange(0.0, 1.0 + bins_step, bins_step)
        self.y_bins_labels = (np.diff(self.y_bins) / 2) + self.y_bins[:-1]
        self.x_bins = np.arange(0.0, 0.5 + bins_step, bins_step) if self.x_data in ['maf_info', 'maf'] \
            else np.arange(0.0, 1.0 + bins_step, bins_step)
        self.x_bins_labels = (np.diff(self.x_bins) / 2) + self.x_bins[:-1]

        # self.quantiles = np.array([0.0, 0.01, 0.1, 0.5, 0.9, 0.99, 1.0])
        self.quantiles = np.array([0.0, 0.01, 0.25, 0.5, 0.75, 0.99, 1.0])

    @property
    def binnedX(self):
        binx = pd.cut(self.dX.values.squeeze(), bins=self.x_bins, labels=self.x_bins_labels)
        return pd.DataFrame(binx, index=self.dX.index, columns=['binned_' + self.x_data])

    @property
    def binnedY(self):
        biny = pd.cut(self.dY.values.squeeze(), bins=self.y_bins, labels=self.y_bins_labels)
        return pd.DataFrame(biny, index=self.dY.index, columns=['binned_' + self.y_data])

    @property
    def quantilY(self):
        df = self.binnedX.join(self.dY)
        pdf = df.groupby(['binned_' + self.x_data]).quantile(self.quantiles)
        pdf = pdf.reset_index()
        pdf.columns = ['binned_' + self.x_data, 'quantiles', self.y_data]
        pdf.astype(float, copy=False)
        return pdf

    @property
    def quantilX(self):
        df = self.binnedY.join(self.dX)
        pdf = df.groupby(['binned_' + self.y_data]).quantile(self.quantiles)
        pdf = pdf.reset_index()
        pdf.columns = ['binned_' + self.y_data, 'quantiles', self.x_data]
        pdf.astype(float, copy=False)
        return pdf

    def rolling_quantilY(self, rollwin: int = 10):
        df = self.dX.join(self.dY).sort_values(by=self.x_data, axis=0)
        rollX = df[self.x_data].rolling(window=rollwin)
        rollY = df[self.y_data].rolling(window=rollwin)

        rollmean = rollX.mean().to_frame()
        dquants = []
        for quant in self.quantiles:
            rollq = rollY.quantile(quant).to_frame()
            rollq['quantiles'] = [quant] * rollq.shape[0]
            rollq = rollmean.join(rollq)
            dquants.append(rollq)

        dfroll = pd.concat(dquants, axis=0).dropna()
        return dfroll

    def binnedX_rolling_quantilY(self, rollwin: int = 10):
        dfroll = self.rolling_quantilY(rollwin=rollwin)
        binx = pd.cut(dfroll[self.x_data].values.squeeze(), bins=self.x_bins, labels=self.x_bins_labels)
        dfrollbin = pd.DataFrame(binx, index=dfroll.index, columns=['binned_' + self.x_data])
        dfrollbin = dfrollbin.join(dfroll[[self.y_data, 'quantiles']])
        df = dfrollbin.groupby(['binned_' + self.x_data, 'quantiles']).mean().reset_index()
        return df



if __name__=='__main__':
    import seaborn as sns
    import matplotlib.pyplot as plt

    true = '/home/camille/1000Genomes/src/VCFPooling/examples/IMP.chr20.snps.gt.vcf.gz'
    imputed_beagle = '/home/camille/1000Genomes/src/VCFPooling/examples/IMP.chr20.pooled.imputed.vcf.gz'
    true = imputed_beagle

    # true = '/home/camille/PoolImpHuman/data/20200722/IMP.chr20.snps.gt.vcf.gz'
    # imputed_beagle = '/home/camille/PoolImpHuman/data/20200722/IMP.chr20.pooled.imputed.vcf.gz'
    # imputed_phaser = '/home/camille/PoolImpHuman/data/20200817/IMP.chr20.pooled.imputed.vcf.gz'

    qbeaglegt = QualityGT(true, imputed_beagle, 0, idx='chrom:pos')

    bgldiff = qbeaglegt.diff()

    # qphasergt = QualityGT(true, imputed_phaser, 0, idx='chrom:pos')
    print(qbeaglegt.trueobj.af_info)

    mafS = qbeaglegt.trueobj.maf_info
    yS_beagle = qbeaglegt.concordance()  # .accuracy

    # pdf_beagle = QuantilesDataFrame(mafS, yS_beagle)
    # pctY_beagle = pdf_beagle.quantilY
    # pctY_beagle['dataset'] = ['beagle'] * pctY_beagle.shape[0]

    dfjoin = mafS.join(yS_beagle).sort_values(by='maf_info', axis=0)
    dfrolled = dfjoin.rolling(window=10).mean().dropna()

    pdf_beagle = QuantilesDataFrame(mafS, yS_beagle)
    pctY_beagle = pdf_beagle.binnedX_rolling_quantilY()  # .quantilY
    pctY_beagle['dataset'] = ['beagle'] * pctY_beagle.shape[0]

    # yS_phaser = qphasergt.concordance()  # .accuracy
    # pdf_phaser = QuantilesDataFrame(mafS, yS_phaser)
    # pctY_phaser = pdf_phaser.quantilY
    # pctY_phaser['dataset'] = ['phaser'] * pctY_phaser.shape[0]
    #
    # pctY_comp = pd.concat([pctY_beagle, pctY_phaser])
    #
    # print(pctY_beagle)

    sns.set(rc={'figure.figsize': (10, 8)})
    sns.set_style('whitegrid')
    dash_styles = [
        (1, 1),
        (3, 1, 1.5, 1),
        (5, 1, 1, 1),
        (5, 1, 2, 1, 2, 1),
        (2, 2, 3, 1.5),
        (1, 2.5, 3, 1.2),
        "",
        (4, 1.5),
    ]

    gY = sns.lineplot(data=pctY_beagle[pctY_beagle.quantiles == 0.5], x='binned_maf_info', y='concordance',
                      palette="deep", linewidth=1)
    gY.fill_between(pctY_beagle[pctY_beagle.quantiles == 1.0]['binned_maf_info'],
                    pctY_beagle[pctY_beagle.quantiles == 0.0]['concordance'],
                    pctY_beagle[pctY_beagle.quantiles == 1.0]['concordance'],
                    color=sns.color_palette('deep')[0],
                    alpha=0.1)
    gY.fill_between(pctY_beagle[pctY_beagle.quantiles == 0.99]['binned_maf_info'],
                    pctY_beagle[pctY_beagle.quantiles == 0.01]['concordance'],
                    pctY_beagle[pctY_beagle.quantiles == 0.99]['concordance'],
                    color=sns.color_palette('deep')[0],
                    alpha=0.25)
    gY.fill_between(pctY_beagle[pctY_beagle.quantiles == 0.9]['binned_maf_info'],
                    pctY_beagle[pctY_beagle.quantiles == 0.1]['concordance'],
                    pctY_beagle[pctY_beagle.quantiles == 0.9]['concordance'],
                    color=sns.color_palette('deep')[0],
                    alpha=0.40)

    plt.show()
