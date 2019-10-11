from cyvcf2 import VCF
import os
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import *
from scripts.poolSNPs import parameters as prm
from scripts.poolSNPs.alleles import alleles_tools as alltls
from scripts.python import results as res

os.chdir(prm.WD)
plots_path = prm.PLOTS_PATH
plt.rcParams["figure.figsize"] = [10, 10]
plt.rcParams['figure.constrained_layout.use'] = True

chk_sz = prm.CHK_SZ
ALL = 'ALL.chr20.snps.gt.chunk{}.vcf.gz'.format(str(chk_sz))

maf_idx, pop_idx = res.make_index(ALL)
data = alltls.vcf2array(VCF(ALL), len(pop_idx))[0]

# MEAN LMAF PER POPULATION
d1 = np.apply_along_axis(alltls.minor_allele, -1, data)
df_maf = pd.DataFrame(data=d1, index=maf_idx)
df_maf.reset_index(level=['maf', 'maf_inter'], drop=False, inplace=True)
df_maf.drop('maf', axis=1, inplace=True)
df_maf.query('maf_inter < 3', inplace=True)
agg_lmaf = df_maf.drop('maf_inter', axis=1, inplace=False)
low_idx = agg_lmaf.index
agg_lmaf = pd.DataFrame(data=agg_lmaf.values, index=low_idx, columns=pop_idx)
agg_lmaf = agg_lmaf.groupby(axis=1, level='Population').mean()
agg_lmaf.boxplot(showfliers=False, showmeans=True)
plt.suptitle('Mean MAF per population for lowest-frequency alleles: True dataset')
plt.savefig(os.path.join(plots_path, 'chunk1000/mean.lmaf.per.pop.png'), dpi='figure')

# MEAN LMAF 1 & 2 PER POPULATION
lev_lmaf = df_maf.set_index('maf_inter', drop=True, append=True, inplace=False)
lev_lmaf = pd.DataFrame(data=lev_lmaf.values, index=lev_lmaf.index, columns=pop_idx)
lev_lmaf = lev_lmaf.groupby(axis=1, level='Population').mean()
lev_lmaf.boxplot(by='maf_inter', showfliers=False, showmeans=True)
plt.suptitle('Mean MAF per population for [0.00-0.01) and [0.01, 0.05) frequencies: True dataset')
plt.savefig(os.path.join(plots_path, 'chunk1000/mean.12maf.per.pop.png'), dpi='figure')

# MEAN HET PER POPULATION
d2 = np.apply_along_axis(alltls.tag_heteroz, -1, data)
df_het = pd.DataFrame(data=d2, index=maf_idx)
df_het.reset_index(level=['maf', 'maf_inter'], drop=False, inplace=True)
df_het.drop('maf', axis=1, inplace=True)
df_het.query('maf_inter < 3', inplace=True)
agg_het = df_maf.drop('maf_inter', axis=1, inplace=False)
low_idx = agg_het.index
agg_het = pd.DataFrame(data=agg_het.values, index=low_idx, columns=pop_idx)
agg_het = agg_het.groupby(axis=1, level='Population').mean()
agg_het.boxplot(showfliers=False, showmeans=True)
plt.suptitle('Mean HET per population for lowest-frequency alleles: True dataset')
plt.savefig(os.path.join(plots_path, 'chunk1000/mean.het.per.pop.png'), dpi='figure')