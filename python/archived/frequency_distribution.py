import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import *
from scipy.optimize import curve_fit, root
from scipy.interpolate import Akima1DInterpolator
from typing import *
import math

from scripts.VCFPooling.poolSNPs import parameters as prm
from scripts.VCFPooling.python.archived.alleles import alleles_tools as alltls
from persotools.files import *


class SigmoidInterpolator(object):
    """

    """
    def __init__(self, file_true_aaf, file_twist_aaf):
        aaf_true = alltls.PandasVCF(file_true_aaf, indextype='chrom:pos').aaf
        aaf_twist = alltls.PandasVCF(file_twist_aaf, indextype='chrom:pos').aaf
        aaf = aaf_true.to_frame().join(aaf_twist, how='inner', rsuffix='_pooled')
        aaf = aaf.sort_values('aaf_pooled', axis=0)
        self.aaf = aaf.drop_duplicates(subset='aaf_pooled', keep='first')

    def get_aaf(self):
        return self.aaf

    @staticmethod
    def sigmoid(x, x0, k, z):
        y = z - (1 / (1 + np.exp(-k * (x - x0))))
        return y

    def fit_sigmoid(self):
        ydata = self.aaf['aaf_pooled'].values.flatten()
        xdata = self.aaf['aaf'].values.flatten()
        popt, pcov = curve_fit(self.sigmoid, xdata, ydata)
        return popt

    def set_sigmoid_params(self):
        popt = self.fit_sigmoid()
        self.__setattr__('params', popt)

    def get_sigmoid_params(self):
        self.set_sigmoid_params()
        return self.__getattribute__('params')

    def call_sigmoid(self, popt: List[float], val: float):
        return self.sigmoid(val, *popt)

    def find_sigmoid_root(self, sol):
        popt = self.fit_sigmoid()
        sigmo_fit = np.vectorize(lambda x: self.sigmoid(x, *popt) - sol)
        r = root(sigmo_fit, [sol])
        return r.x[0]

    def vectorize_sigmoid_root(self):
        nb_pts = 1000
        x_sol = np.arange(0, 1.0 + (1/nb_pts), 1/nb_pts)
        y_sol = np.clip([self.find_sigmoid_root(x) for x in x_sol],
                        0.0,
                        1.0)
        return x_sol, y_sol

    def sigmoid_derivative(self):
        self.set_sigmoid_params()
        sigvect = np.vectorize(lambda val: self.call_sigmoid(self.get_sigmoid_params(), val))
        nb_pts = 1000
        dx = 1/nb_pts
        x_val = np.arange(0, 1.0 + dx, dx)
        y_val = sigvect(x_val)
        yy_val = np.gradient(y_val, dx)
        return x_val, yy_val

    def interpolate_sigmoid(self):
        xi, yi = self.vectorize_sigmoid_root()
        akima = Akima1DInterpolator(xi, yi)
        return akima

    def interpolate_derivative(self):
        xi, yyi = self.sigmoid_derivative()
        akima = Akima1DInterpolator(xi, yyi)
        return akima

    def call_reciprocal_sigmoid(self, val: float):
        interpolator = self.interpolate_sigmoid()
        return interpolator.__call__(np.array(val))

    @staticmethod
    def call_sigmoid_derivative(interpolator: Any, val: float):
        return interpolator.__call__(np.array(val))


def write_sigmoid_csv(nb_points=50):
    os.chdir(prm.PATH_GT_FILES)

    sig = SigmoidInterpolator('ALL.chr20.snps.gt.chunk10000.vcf.gz',
                              'ALL.chr20.pooled.snps.gt.chunk10000.vcf.gz')
    sig.set_sigmoid_params()
    aaf = sig.get_aaf()

    x_, y_ = sig.vectorize_sigmoid_root()
    params = sig.get_sigmoid_params()

    reciprocvect = np.vectorize(lambda x: sig.call_reciprocal_sigmoid(x))
    sigvect = np.vectorize(lambda x: sig.call_sigmoid(params, x))
    x_pts = np.arange(0.0, 1.0, 1.0/nb_points)
    y_pts = sigvect(x_pts)
    y_1_pts = reciprocvect(x_pts)
    df_sig = pd.DataFrame(zip(x_pts, y_pts, y_1_pts),
                          columns=['x', 'y_sigmoid', 'y_reciprocal'])
    df_sig.to_csv(os.path.join(prm.DATA_PATH,
                               'sigmoid_coordinates_{}points.csv'.format(nb_points)),
                  sep=',',
                  index=False)


def gl_distribution():
    rr = np.vectorize(lambda x: (1 - x)**2)
    ra = np.vectorize(lambda x: 2 * (1 - x) * x)
    aa = np.vectorize(lambda x: x ** 2)
    return rr, ra, aa


def convolve_gl(nb_points=50):
    f = np.arange(0.0, 1.0, 1.0/nb_points)
    sig = pd.read_csv(os.path.join(prm.DATA_PATH,
                                   'sigmoid_coordinates_{}points.csv'.format(nb_points)),
                      sep=',').loc[:, 'y_sigmoid']
    ones = np.ones_like(f)
    sig_ones = np.subtract(sig, ones)
    rr, ra, aa = gl_distribution()
    compo_rr = rr(sig)
    compo_ra = ra(sig)
    compo_aa = aa(sig)
    convo_rr = np.convolve(sig, f, mode='same')  # alltls.normalize(convo_rr)
    plt.plot(f, rr(f), 'b--', label='RR')
    plt.plot(f, compo_rr, 'b-', label='composed RR')
    plt.plot(f, ra(f), 'g--', label='RA')
    plt.plot(f, compo_ra, 'g-', label='composed RA')
    plt.plot(f, aa(f), 'r--', label='AA')
    plt.plot(f, compo_aa, 'r-', label='composed AA')
    plt.legend()
    plt.show()


def plot_sigmoid_aaf(plot=False):
    os.chdir(prm.PATH_GT_FILES)

    sig = SigmoidInterpolator('ALL.chr20.snps.gt.chunk10000.vcf.gz',
                              'ALL.chr20.pooled.snps.gt.chunk10000.vcf.gz')
    sig.set_sigmoid_params()
    aaf = sig.get_aaf()

    x_, y_ = sig.vectorize_sigmoid_root()
    params = sig.get_sigmoid_params()
    interp = sig.interpolate_derivative()
    g_hetero = lambda x: 2 * x * (1-x)
    g_hetero2 = lambda x: 4 * x * (1-x)

    x_point = 0.48
    y_point = sig.call_sigmoid(params, x_point)
    twist_ratio = x_point/y_point
    print('f={}, f_twisted={}, f/f_twist={}'.format(x_point, y_point, twist_ratio))
    gl = math.log10(g_hetero(x_point))
    # x_pow4 = x_point * (twist_ratio ** 4)
    x_pow4 = gl + math.log10(twist_ratio ** 4)
    print(x_pow4)
    gl_pow4 = math.log10(g_hetero(x_pow4))
    gl2 = math.log10(g_hetero2(x_point))
    print('gl={}, gl_pow4={}, gl2={}'.format(gl, gl_pow4, gl2))

    if plot:
        sigvect = np.vectorize(lambda x: sig.call_sigmoid(params, x))
        x_pts = np.arange(0.0, 1.0, 0.01)
        y_pts = sigvect(x_pts)
        y_x = np.fabs(np.divide(np.subtract(y_pts, x_pts), x_pts))
        yx = np.divide(y_pts, x_pts)
        deriv = np.vectorize(lambda x: sig.call_sigmoid_derivative(interp, x))
        yy_pts = deriv(x_pts)

        fig, ax = plt.subplots()
        # ax.plot(x_, y_, linestyle='-', c='b', label='fit')
        # ax.scatter(aaf.values[:, 1], aaf.values[:, 0], marker='o', s=2, c='g', label='true')
        ax.plot(x_pts, y_pts, linestyle='-', c='r', label='points')
        ax.plot(x_pts, x_pts, linestyle='-', c='k', label='identity')
        # ax.plot(x_pts, y_x, linestyle='--', c='k', label=' difference aat-aaf')
        ax.plot(x_pts, yx, linestyle='--', c='k', label='ratio')
        ax.plot(x_pts, yy_pts, linestyle='--', c='r', label='derivative')
        ax.legend()
        plt.show()


def expected_decodability(nb_samples: int = 16, fqc: np.ndarray = np.linspace(0.0, 1.0, 100), d: int = 1):
    """

    :param nb_samples: number of items in the pool
    :param fqc: array of ALT allele frequencies
    :param d: decoding robustness
    :return:
    """
    e_pos = nb_samples * fqc
    plt.plot(fqc, e_pos, color='b', label='expected number of ALT carriers per pool')
    plt.axhline(y=d, xmin=0.0, xmax=1.0, color='k', linestyle='--', label='decoding robustness')
    plt.xlabel('ALT allele frequency')
    plt.show()


if __name__ == '__main__':
    # write_sigmoid_csv()
    # convolve_gl()
    #plot_sigmoid_aaf()
    expected_decodability()
