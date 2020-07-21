from scripts.VCFPooling.poolSNPs.poolvcf import *
from scripts.VCFPooling.python.archived.pool import *
import numpy as np
import timeit


def time_process_file(data: VCF, groups: list, simul: str, vcf_out: str) -> None:
    """
    Computes and rewrites genotypes of all individuals for all samples from input files
    :param data: cyvcf2 object reader pointing on a VCF-file
    :param groups: samples identifiers split in pools
    :param f: integer, index of the file to process in the list
    :param fileout: VCF-files with simulated pooled or randomly missing genotypes
    """
    print('Simulation type: ', simul)
    print('File out: ', os.path.join(os.getcwd(), vcf_out))  # prm.PATH_OUT[simul]
    if prm.GTGL == 'GL' and prm.unknown_gl == 'adaptive':
        dic_header = {'ID': 'GL',
                      'Number': 'G',
                      'Type': 'Float',
                      'Description': 'three log10-scaled likelihoods for RR,RA,AA genotypes'}
        data.add_format_to_header(dic_header)
        whead = Writer(vcf_out, data)
        whead.write_header()
        whead.close()
        w = open(vcf_out, 'ab')
        # Load adaptive GL values for missing data
        df = pd.read_csv(os.path.join(prm.WD, 'adaptive_gls.csv'),
                         header=None,
                         names=['rowsrr', 'rowsra', 'rowsaa', 'colsrr', 'colsra', 'colsaa',
                                'n', 'm',
                                'rr', 'ra', 'aa']
                         )
        df2dict = dict(((int(rwrr), int(rwra), int(rwaa), int(clrr), int(clra), int(claa),
                         int(n), int(m)),
                        [rr, ra, aa]) for rwrr, rwra, rwaa, clrr, clra, claa,
                                          n, m,
                                          rr, ra, aa in df.itertuples(index=False, name=None))

    else:  # prm.GTGL == 'GT' or fixed GL
        w = Writer(vcf_out, data)
        w.set_threads(4)
        df2dict = None

    tm = timeit.default_timer()
    # for n, variant in enumerate(data('20:59973567-59973568')):
    for n, variant in enumerate(data):
        process_line(groups, simul, w, variant, df2dict, write=True)
        if n % 1000 == 0:
            print('{} variants processed in {:06.2f} sec'.format(n+1, timeit.default_timer()-tm).ljust(80, '.'))
        # if n+1 == 1000:
        #     break
    w.close()

    # GL converted from GT, missing GLs will be filled with [0.33, 0.33, 0.33]
    if prm.GTGL == 'GL' and prm.unknown_gl != 'adaptive':
        alltls.file_likelihood_converter(os.path.join(prm.PATH_GT_FILES,
                                                      vcf_out.replace('.gl', '.gt')) + '.gz',
                                         vcf_out)

datadir = '/home/camille/PoolImpHuman/data/20200603'
os.chdir(datadir)


### Last version of implementation for simulating pooling on a VCF file
design = Design()
dm = design.matrix

path_keys = '/home/camille/PoolImpHuman/data/main'
dict_gl = load_lookup_dict(path_keys)

poolf = VariantFilePooler(dm,
                          os.path.join(datadir, "TEST.chr20.snps.gt.vcf.gz"),
                          os.path.join(datadir, "TEST.chr20.snps.pooled.gp.vcf"),
                          dict_gl,
                          'GP')

tstart = timeit.default_timer()
poolf.write()
tstop = timeit.default_timer()
print('Time for pooling {} variants = {} sec'.format(poolf.n_variants, tstop-tstart))
"""
1 variants processed in 000.17 sec..............................................
"""

### Older version
vcf = VCF(os.path.join(os.getcwd(), "TEST.chr20.snps.gt.vcf.gz"))
samples = poolf.vcf_in.header.samples
n_splits = len(samples) // poolf.design.shape[1]
splits = np.array(samples).reshape((n_splits, poolf.design.shape[1]))

# run_pool.pt output:
"""
Number of 16-sized pools to create, number of single samples remaining:  156 8
/home/camille/1000Genomes/data/gt/ALL.chr20.snps.gt.chunk1000.unordered.samples.vcf.gz:
 File created? -> True
/home/camille/1000Genomes/data/gt/ALL.chr20.snps.gt.chunk1000.vcf.gz:
 File indexed? -> True
/home/camille/1000Genomes/data/gt/ALL.chr20.snps.gt.chunk1000.vcf.gz:
 File created? -> True
/home/camille/1000Genomes/data/gt/ALL.chr20.snps.gt.chunk1000.vcf.gz:
 File indexed? -> True
ALL.chr20.snps.gt.chunk1000.unordered.samples.vcf.gz.csi: the file does not exists
File to write to -->  ALL.chr20.pooled.snps.gl.chunk.chunk1000.vcf
Simulation type:  pooled
File out:  /home/camille/1000Genomes/data/gl/ALL.chr20.pooled.snps.gl.chunk.chunk1000.vcf
1 variants processed in 000.94 sec..............................................
"""
