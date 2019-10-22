import numpy as np
import matplotlib.pyplot as plt
from cyvcf2 import VCF

from scripts.VCFPooling.poolSNPs.alleles import alleles_tools as alltls
from scripts.VCFPooling.poolSNPs import parameters as prm
from persotools.files import *

"""
RGB proportional to GL(RR|RA|AA)
"""

def gl_colors(gen):
    """

    :param gen: genotype (log-)likelihoods. Array-like.
    :return: triplet coding rgb
    """
    tenpow = np.vectorize(lambda x: 0.0 if x == -5.0 else pow(10, x))
    rgb = tenpow(gen)
    return rgb


def gt_colors(gen):
    """

    :param gen: raw genotypes with phase. Array-like.
    :return: triplet coding rgb
    """
    rgb = np.apply_along_axis(alltls.map_gt_gl, axis=-1, arr=gen)
    return rgb


def gl2gt(gen):
    """
    Returns the most likely genotype.
    :param gen: triplet for genotype likelihoods
    :return:
    """
    if gen[0] == 1.0:
        gt = '0/0'
    elif gen[2] == 1.0:
        gt = '1/1'
    elif gen[1] == 1.0:
        gt = '1/0'
    elif round(gen[1], 1) == 0.5 and round(gen[2], 1) == 0.5:
        gt = '1/.'
    elif round(gen[0], 1) == 0.5 and round(gen[1], 1) == 0.5:
        gt = '0/.'
    else:
        gt = './.'
    return gt


def read_gl_vcf(file_path):
    """

    :param file_path: file to reade
    :return:
    """
    import subprocess
    subprocess.run('bcftools view -H {} > {}'.format(file_path, 'temp.gl.txt'),
                   shell=True,
                   cwd=os.getcwd())
    with open('temp.gl.txt', mode='r') as vcf:
        content = vcf.readlines()

    for line in content:
        s = line
        if s.startswith('##'):
            continue
        elif s.startswith('#CHROM'):
            header = s
        else:
            s = s[:-1] # trim '\n'
            var = s.split('\t')
            pos = var[1:2]
            gl = var[9:]
            if len(pos) > 0:
                yield pos[0], [loggl.split(',') for loggl in gl]
    delete_file('temp.gl.txt')


def plot_pools(vector, gtgl, step, snp, sz, af_info=None, aaf=None):
    """

    :param marker: log-GL of samples, str format. NumPy array.
    :return:
    """
    if af_info is None:
        af_info = '.'
    if aaf is None:
        aaf = '.'

    plt.rcParams['axes.titlesize'] = 6
    plots_sz = (3, 5)
    fig, axes = plt.subplots(plots_sz[0], plots_sz[1],
                             figsize=(4, 4))
    fig.suptitle('SNP 20:{}; AF_INFO = {}; AAF = {}'.format(snp, af_info, aaf),
                 fontsize=6)

    if gtgl == 'gl':
        pool = np.apply_along_axis(gl_colors, axis=-1, arr=vector)
    if gtgl == 'gt':
        pool = np.apply_along_axis(gt_colors, axis=-1, arr=vector)
    tx = np.apply_along_axis(gl2gt, axis=-1, arr=pool)
    pool = np.multiply(pool, np.broadcast_to([0.5, 0.8, 0.5], pool.shape))  # modify RGB colors
    k = 0
    for i in range(1, plots_sz[0] + 1):
        for j in range(1, plots_sz[1] + 1):
            k += 1
            axes[i - 1, j - 1].imshow(pool[k-1, :, :].reshape((4, 4, 3)),
                                      cmap='plasma')
            axes[i - 1, j - 1].set_title('Pool #{}'.format(k), pad=2)
            axes[i - 1, j - 1].set_xticks(np.arange(4) + 0.5, minor=True)
            axes[i - 1, j - 1].set_yticks(np.arange(4) + 0.5, minor=True)
            axes[i - 1, j - 1].tick_params(which="both",
                                           axis='both',
                                           bottom=False,
                                           left=False,
                                           labelsize=4,
                                           length=1,
                                           pad=0.3)
            axes[i - 1, j - 1].grid(which='minor', axis='both',
                                    color="w", linestyle='-', linewidth=0.25)
            # remnove borders
            axes[i - 1, j - 1].spines['top'].set_visible(False)
            axes[i - 1, j - 1].spines['right'].set_visible(False)
            axes[i - 1, j - 1].spines['bottom'].set_visible(False)
            axes[i - 1, j - 1].spines['left'].set_visible(False)
            tx_i_j = tx[k - 1, :].reshape((4, 4))
            for m in range(4):
                for n in range(4):
                    axes[i - 1, j - 1].text(n, m, tx_i_j[m, n], ha="center", va="center", color="w", fontsize=4)
    plt.autoscale()
    fig.tight_layout()
    plt.savefig(prm.PLOTS_PATH + '/pools_patterns.{}.chunk{}.snp{}.jpg'.format(step, sz, snp),
                dpi=500)


def before_after_pooling(snp, sz):
    """
    Combine pools overview in a PDF document
    :param snp:
    :param sz:
    :return:
    """
    import img2pdf
    imglist = [prm.PLOTS_PATH + '/pools_patterns.0.chunk{}.snp{}.jpg'.format(sz, snp),
               prm.PLOTS_PATH + '/pools_patterns.1.chunk{}.snp{}.jpg'.format(sz, snp),
               prm.PLOTS_PATH + '/pools_patterns.2.chunk{}.snp{}.jpg'.format(sz, snp)]
    with open(prm.PLOTS_PATH + '/pools_patterns.chunk{}.snp{}.pdf'.format(sz, snp), "wb") as f:
        f.write(img2pdf.convert([i for i in imglist]))


# TODO: which frequency does Beagle use?
def plot_af_vs_aaf(file_all, file_imp, file_ref, file_pool_imp, file_pool_ref, sz):
    """
    Show variations/differences between:
        - AAF estimated on the whole dataset (hardcoded in the VCF file);
        - AAF computed on the datasets itself (ALL, IMP, REF)
    :return:
    """
    allaafs = alltls.get_aaf(file_all,
                             id='chrom:pos')
    impaafs = alltls.get_aaf(file_imp,
                             id='chrom:pos')
    refaafs = alltls.get_aaf(file_ref,
                             id='chrom:pos')
    poolimp = alltls.get_aaf(file_pool_imp,
                             id='chrom:pos')
    poolref = alltls.get_aaf(file_pool_ref,
                             id='chrom:pos')
    compaafs = allaafs.loc[:, ['af_info', 'aaf']].join(impaafs.loc[:, ['af_info', 'aaf']], how='inner', rsuffix='_imp')
    compaafs = compaafs.join(refaafs.loc[:, ['af_info', 'aaf']], how='inner', rsuffix='_ref')
    compaafs = compaafs.join(poolimp.loc[:, ['af_info', 'aaf']], how='inner', rsuffix='_pooled_imp')
    compaafs = compaafs.join(poolref.loc[:, ['af_info', 'aaf']], how='inner', rsuffix='_pooled_ref')

    compaafs.sort_values(['af_info', 'aaf_imp'], axis=0, inplace=False)
    print(compaafs)

    fig, axis = plt.subplots()
    compaafs.plot.scatter(x='af_info',
                          y='af_info',
                          s=3,
                          ax=axis,
                          label='AAF from VCF AF INFO-field',
                          marker='o',
                          color='tab:blue')
    compaafs.plot.scatter(x='af_info',
                          y='aaf',
                          s=3,
                          ax=axis,
                          label='AAF from VCF variants for all samples',
                          marker='o',
                          color='tab:green')
    compaafs.plot.scatter(x='af_info',
                          y='aaf_imp',
                          s=3,
                          ax=axis,
                          label='AAF from VCF variants for imputation samples without pooling-decoding',
                          marker='o',
                          color='tab:orange')
    compaafs.plot.scatter(x='af_info',
                          y='aaf_ref',
                          s=3,
                          ax=axis,
                          label='AAF from VCF variants for reference panel without pooling-decoding',
                          marker='o',
                          color='tab:purple')
    compaafs.plot.scatter(x='af_info',
                          y='aaf_pooled_imp',
                          s=3,
                          ax=axis,
                          label='AAF from VCF variants for imputation samples after pooling-decoding',
                          marker='o',
                          color='tab:red')
    compaafs.plot.scatter(x='af_info',
                          y='aaf_pooled_ref',
                          s=3,
                          ax=axis,
                          label='AAF from VCF variants for reference panel after pooling-decoding',
                          marker='o',
                          color='tab:pink')
    axis.legend(loc='center', fontsize=6, bbox_to_anchor=(0.5, -0.3))
    plt.xlabel('AAF from VCF AF INFO-field')
    plt.ylabel('Alternate allele frequency')
    fig.tight_layout()
    plt .savefig(prm.PLOTS_PATH + '/pooled_af_vs_aaf.chunk{}.jpg'.format(sz),
                dpi=500)


if __name__ == '__main__':
    ### DATA
    if prm.GTGL == 'GT':
        cd = prm.PATH_GT_FILES
    if prm.GTGL == 'GL':
        if prm.unknown_gl != 'adaptative':
            cd = os.path.join(prm.WD, 'gl', 'gl_' + '_'.join(np.around(prm.unknown_gl, decimals=2).astype(str)))
        else:
            # cd = os.path.join(prm.WD, 'gl')
            cd = os.path.join(prm.WD, 'gl', 'gl_adaptative')

    print('Load parameters in {}'.format(cd).ljust(80, '.'))
    sorting = True  # sort data sets by AAF and population values
    os.chdir(cd)
    chk_sz = prm.CHK_SZ

    params = ('gl', chk_sz, 'chrom:pos')
    devol = []

    # Load AAFs
    aafs = alltls.get_aaf(os.path.join(prm.PATH_GT_FILES,
                                       'IMP.chr20.pooled.snps.gt.chunk{}.vcf.gz'.format(chk_sz)),
                          id='chrom:pos')
    aafs = aafs.loc[:, 'af_info'].to_frame()
    ALL = os.path.join(prm.PATH_GT_FILES, 'ALL.chr20.snps.gt.chunk{}.vcf.gz'.format(chk_sz))
    REF = os.path.join(prm.PATH_GT_FILES, 'REF.chr20.snps.gt.chunk{}.vcf.gz'.format(chk_sz))
    IMP = os.path.join(prm.PATH_GT_FILES, 'IMP.chr20.snps.gt.chunk{}.vcf.gz'.format(chk_sz))
    B1 = 'ALL.chr20.pooled.snps.gl.chunk{}.vcf.gz'.format(chk_sz)
    B1_IMP = os.path.join(prm.PATH_GT_FILES, 'IMP.chr20.pooled.snps.gt.chunk{}.vcf.gz'.format(chk_sz))
    B1_REF = os.path.join(prm.PATH_GT_FILES, 'REF.chr20.pooled.snps.gt.chunk{}.vcf.gz'.format(chk_sz))
    OUT = os.path.join(cd,
                       # 'single_samples_merged',
                       # 'all_snps_all_samples',
                       'phaser',
                       'IMP.chr20.pooled.imputed.gt.chunk{}.vcf.gz'.format(chk_sz))
    # plot_af_vs_aaf(ALL, IMP, REF, B1_IMP, B1_REF, chk_sz)

    candidates = alltls.extract_variant_onfreq(IMP, [0.499, 0.501])
    print(candidates)

    # from leave_one_out vs one-batch experiment
    outliers = ['20:7696268', '20:7145990', '20:8167370', '20:7579308', '20:7495442', '20:8964609', '20:8633433',
                '20:8463605', '20:9765652', '20:7670926', '20:7673476', '20:7445132', '20:7677668', '20:7175142']

    varpos = candidates['id'].str.strip(' ').str.slice(3).values
    # varpos = list([var[3:].strip(' ') for var in outliers])
    # varpos = ['59973567']  # 10000 markers
    # varpos = ['25344231', '42715255', '29827638', '58905133', '38406028', '12617340']  # 1000 markers

    if True:
        for v in varpos:
            print('\nSet params: {} for variant 20:{}'.format(params, v))
            gtgl, sz, idt = params
            varid = ''.join(['20:', v])
            print('AAF: ', aafs.loc[varid])  # id.ljust(11, ' ')
            gt_iterator = VCF(ALL)
            gl_iterator = read_gl_vcf(B1)
            out_iterator = VCF(OUT)
            for g in gt_iterator:
                if str(g.POS) == v:
                    marker0 = np.array(g.genotypes).reshape((156, 16, 3)).astype(float)
                    aaf0 = g.aaf
            for g in gl_iterator:
                if g[0] == v:
                    marker1 = np.array(g[-1]).reshape((156, 16, 3)).astype(float)
            for g in out_iterator:
                if str(g.POS) == v:
                    marker2 = np.array(g.genotypes).reshape((15, 16, 3)).astype(float)
                    aaf2 = g.aaf

            plot_pools(marker0, 'gt', 0, v, chk_sz, af_info=aafs.loc['20:' + v, 'af_info'], aaf=aaf0)
            plot_pools(marker1, 'gl', 1, v, chk_sz, af_info=aafs.loc['20:' + v, 'af_info'])
            plot_pools(marker2, 'gt', 2, v, chk_sz, af_info=aafs.loc['20:' + v, 'af_info'], aaf=aaf2)

            before_after_pooling(v, chk_sz)


