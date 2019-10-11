import subprocess
import numpy as np

from scripts.poolSNPs import parameters as prm
from persotools.files import *

"""
Bash commands for bcftools manipulations written as Python-functions.
Increase readability.
"""


def bgzip(f_vcf: str, f_gz: str, wd: str) -> None:
    """
    Bgzip a vcf file into a vcf.gz and checks if creation succeeded.
    :param f_vcf: input vcf file name
    :param f_gz: output vcf.gz file name
    :param wd: path to working directory
    :return: None
    """
    cmd = ' '.join(['bcftools',
                    'view -Oz -o',
                    f_gz,
                    f_vcf
                    ])
    subprocess.run(cmd, shell=True, cwd=wd)
    print('{}:\r\n File created? -> {}'.format(os.path.join(wd, f_gz),
                                               check_file_creation(wd, f_gz)))


def sort(f_gz: str, wd: str) -> None:
    """
    Sorts a VCF-file per increasing marher position on chromosom.
    :param f_gz: input file name (vcf or vcf.gz)
    :param wd: path to working directory
    :return: None
    """
    cmd = ' '.join(['bcftools',
                    'sort -Oz -o',
                    f_gz,
                    f_gz
                    ])
    subprocess.run(cmd, shell=True, cwd=wd)
    print('{}:\r\n File sorted? -> {}'.format(os.path.join(wd, f_gz),
                                              check_file_creation(wd, f_gz)))


def index(f_gz: str, wd: str) -> None:
    """
    Index a VCF-file and replace old .csi file if it exixts.
    :param f_gz: input file name (vcf or vcf.gz)
    :param wd: path to working directory
    :return: None
    """
    cmd = ' '.join(['bcftools',
                    'index -f',
                    f_gz
                    ])
    subprocess.run(cmd, shell=True, cwd=wd)
    print('{}:\r\n File indexed? -> {}'.format(os.path.join(wd, f_gz),
                                               check_file_creation(wd, f_gz + '.csi')))


def sampling(f_gz: str, f_out: str, f_samp: str, wd: str) -> None:
    """
    Creates a VCF-file (vcf.gz) for a chosen set of samples.
    Samples should be stored in a text-file, one sample name per line.
    Checks if sampled file was created.
    :param f_gz: input VCF-file name
    :param f_out: output VCF-file name (vcf.gz)
    :param f_samp: samples file name (text file)
    :param wd: path to working directory
    :return: None
    """
    cmd = ' '.join(['bcftools',
                    'view -Oz -o',
                    f_out,
                    '-S {}'.format(f_samp),
                    f_gz
                    ])
    subprocess.run(cmd, shell=True, cwd=wd)
    print('{}:\r\n File created? -> {}'.format(os.path.join(wd, f_out),
                                               check_file_creation(wd, f_out)))


def stratified_aaf_sampling(f_gz: str, wd: str, binning: bool = False) -> None:
    """
    Samples markers from a VCF-file by selecting an equal number of
    markers for each 0.1-MAF interval. Total number of markers depends
    on the global chunk size parameter.
    :param f_gz: input VCF-file with all markers
    :param wd: path to working directory
    :param binning: boolean for splitting the original data set into AAF-binned layers
    :return: None
    """
    os.chdir(wd)
    print(wd)
    subprocess.run(' '.join(['bcftools',
                             'view -h -Ov -o',
                             'headers.ALL.chr20.snps.gt.chunk{}.strat.vcf'.format(prm.CHK_SZ),
                             os.path.join(prm.DATA_PATH, 'gt', 'ALL.chr20.snps.gt.vcf.gz')
                             ]),
                   shell=True, cwd=wd)
    bins = np.arange(0.0, 1.0, 0.1)
    for b in bins:
        print('interval for AAF: [{}, {})'.format(b, b+0.1))
        cmd1 = ' '.join(['bcftools',
                         'view -Oz -o',
                         f_gz.replace('.vcf.gz', '.maf_{}_{}.vcf.gz'.format(b, b+0.1)),
                         '-q {} -Q {}'.format(b, b+0.1),
                         f_gz
                         ])

        tmp = ' '.join(['cat ./bins_for_chunk10000/chunk{}.vcf'.format(b),
                        '| sort -R | head -{}'.format(prm.CHK_SZ // len(bins)),
                        '> chunk{}.vcf'.format(b)
                        ])  # for generating strat_chunk1000 from strat_chunk10000

        cmd2 = ' '.join(['bcftools view -H',
                         f_gz.replace('.vcf.gz', '.maf_{}_{}.vcf.gz'.format(b, b+0.1)),
                         '| sort -R | head -{}'.format(prm.CHK_SZ // len(bins)),
                         '> chunk{}.vcf'.format(b)
                         ])

        if binning:
            subprocess.run(cmd1, shell=True, cwd=wd)
        subprocess.run(tmp, shell=True, cwd=wd)
        # subprocess.run(cmd2, shell=True, cwd=wd)

    subprocess.run(' '.join(['cat headers.ALL.chr20.snps.gt.chunk{}.strat.vcf '.format(prm.CHK_SZ),
                             ' '.join(['chunk{}.vcf'.format(i) for i in bins]),
                             '> TMP.chr20.snps.gt.strat.vcf']),
                   shell=True,
                   cwd=wd)

    delete_file('ALL.chr20.snps.gt.chunk{}.vcf.gz.csi'.format(prm.CHK_SZ))

    bgzip('TMP.chr20.snps.gt.strat.vcf'.format(prm.CHK_SZ),
          'ALL.chr20.snps.gt.chunk{}.vcf.gz'.format(prm.CHK_SZ),
          wd)
    sort('ALL.chr20.snps.gt.chunk{}.vcf.gz'.format(prm.CHK_SZ),
         wd)
    index('ALL.chr20.snps.gt.chunk{}.vcf.gz'.format(prm.CHK_SZ),
          wd)
    delete_file('headers.ALL.chr20.snps.gt.chunk{}.strat.vcf'.format(prm.CHK_SZ))
    print('Stratification ended\n\n'.ljust(80, '.'))
    # delete_file('TMP.chr20.snps.gt.strat.vcf')
    # for i in bins:
    #     delete_file('chunk{}.vcf'.format(i))


def rename_samples(file_in: str, file_out:str, wd: str, suffix:str) -> None:
    """
    Rename samples of a file by adding suffix to the old names.
    :param file_in:
    :param file_out:
    :param suffix:
    :return:
    """
    cmd0 = ' '.join(['bcftools query -l',
                     '-l',
                     file_in,
                     '> tmp.samples.get_names.txt'])
    subprocess.run(cmd0, shell=True, cwd=wd)

    with open(os.path.join(wd, 'tmp.samples.get_names.txt'), 'r') as get_names:
        names = get_names.readlines()

    with open(os.path.join(wd, 'tmp.samples.set_names.txt'), 'w') as set_names:
        lines = ['{} {}'.format(n.strip('\n\r'), n.strip('\n\r') + suffix + '\n') for n in names]
        set_names.writelines(lines)

    cmd1 = ' '.join(['bcftools reheader',
                    '-s tmp.samples.set_names.txt',
                     '-o',
                     file_out,
                     file_in])
    subprocess.run(cmd1, shell=True, cwd=wd)

    # delete_file(os.path.join(wd, 'tmp.samples.get_names.txt'))
    # delete_file(os.path.join(wd, 'tmp.samples.set_names.txt'))


if __name__ == '__main__':
    pass
    # os.chdir(prm.WD + '/gt')
    # print(os.getcwd())
    # mkdir(os.path.join(prm.DATA_PATH,
    #                    'gt',
    #                    'stratified'))
    # stratified_aaf_sampling(prm.PATH_IN,
    #                         os.path.join(prm.DATA_PATH,
    #                                      'gt',
    #                                      'stratified'))
