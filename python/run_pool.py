import os
from cyvcf2 import VCF
from scripts.VCFPooling.poolSNPs import parameters as prm
from scripts.VCFPooling.poolSNPs import pool

nb_cores = os.cpu_count()

# TODO: exchange names source - path_in
if prm.GTGL == 'GT':
    CD = prm.PATH_GT_FILES
if prm.GTGL == 'GL':
    CD = prm.PATH_GL_FILES
os.chdir(CD)
CHK_SZ = prm.CHK_SZ
SRCFILE = prm.SRCFILE
PATH_OUT = prm.PATH_OUT
POOL = prm.POOL
source = prm.CHKFILE
# source = 'ALL.chr20.snps.gt.chunk10000.strat.vcf.gz'
subset = prm.SUBSET

groups = pool.init_chunk(CD, SRCFILE, chunk=False, strat=True)

data = VCF(os.path.join(prm.WD, 'gt', source), threads=nb_cores)

SAMPLES = data.samples

pool.run(groups, range(2))
# run(groups, [1])

if subset:
    pool.subset_chunked_vcf(CD, source, PATH_OUT, CHK_SZ, prm.SUBCHUNK)