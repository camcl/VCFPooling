import os
from cyvcf2 import VCF
from src.VCFPooling.poolSNPs import parameters as prm
from src.VCFPooling.python.archived import pool

nb_cores = os.cpu_count()

if prm.GTGL == 'GT':
    CD = prm.PATH_GT_FILES
if prm.GTGL == 'GL':
    CD = prm.PATH_GL_FILES
os.chdir(CD)
CHK_SZ = prm.CHK_SZ
SRCFILE = prm.SRCFILE
PATH_OUT = prm.PATH_OUT
POOL = prm.POOL
chkfile = prm.CHKFILE
subset = prm.SUBSET

groups = pool.init_chunk(CD, SRCFILE, chunk=False, strat=True)

data = VCF(os.path.join(prm.WD, 'gt', chkfile), threads=nb_cores)

SAMPLES = data.samples

for sim in prm.PATH_OUT.keys():
    pool.run(groups, sim)

if subset:
    pool.subset_chunked_vcf(CD, chkfile, PATH_OUT, CHK_SZ, prm.SUBCHUNK)
