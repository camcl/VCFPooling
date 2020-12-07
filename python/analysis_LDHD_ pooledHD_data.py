import os, sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

from typing import *

rootdir = os.path.dirname(os.path.dirname(os.getcwd()))
sys.path.insert(0, rootdir)

from VCFPooling.poolSNPs.metrics import quality as qual
from VCFPooling.poolSNPs import dataframe as vcfdf
from VCFPooling.persotools.files import *

'''
Are all the correctly decoded markers/fully assayed LD markers remaining unchanged during imputation?
What is the error on undecoded/unasseyed HD markers after imputation?
Is there a correlation between the markers linkage and what data/how much is missing?
Look at a single sample HG01063 (Phaser works on single samples anyways)
'''


# Data parameters and plotting features

sns.set(rc={'figure.figsize': (10, 8)})
sns.set_style('whitegrid')


# Configure data/plots paths

outdir = '/home/camille/PoolImpHuman/results/20201029'
if not os.path.exists(outdir):
    os.mkdir(outdir)

truegt = '/home/camille/PoolImpHuman/data/20201029/sHG01063.IMP.chr20.snps.gt.vcf.gz'
truegl = '/home/camille/PoolImpHuman/data/20201029/sHG01063.IMP.chr20.snps.gl.vcf.gz'
# pooled can also be the file with full LD and missing HD
pooledgt = '/home/camille/PoolImpHuman/data/20201029/sHG01063.IMP.chr20.missingHD.fullLD.snps.gt.vcf.gz'
pooledgl = None
# imputation with Beagle or with Phaser
imputed_gtgp1 = '/home/camille/PoolImpHuman/data/20201029/sHG01063.IMP.chr20.pooled.imputed.vcf.gz'
imputed_gtgp2 = None


# Build Quality and DataFrame objects for analysis

quality1gt = qual.QualityGT(truegt, imputed_gtgp1, 0, idx='chrom:pos')
quality1gl = qual.QualityGL(truegl, imputed_gtgp1, 0, idx='chrom:pos')

if imputed_gtgp2 is not None:
    quality2gt = qual.QualityGT(truegt, imputed_gtgp2, 0, idx='chrom:pos')
    quality2gl = qual.QualityGL(truegl, imputed_gtgp2, 0, idx='chrom:pos')

if pooledgt is not None:
    dfpooled = vcfdf.PandasMixedVCF(pooledgt, format='GT', indextype='chrom:pos')

if pooledgl is not None:
    dfpooled = vcfdf.PandasMixedVCF(pooledgl, format='GL', indextype='chrom:pos')


# Coordinates of HDonly and LDonly variants chr20 x Illumina (35,682 and 17,015 variants)

ld_vars = pd.read_csv('/home/camille/PoolImpHuman/data/omniexpress-isec-LDHD/LDonly.coords',
                         sep='\t',
                         header=None,
                         names=['CHROM', 'POS'])
ld_vars['variants'] = pd.Series(['20:{}'.format(pos) for pos in ld_vars['POS']], dtype=str)
ld_vars.set_index('variants', inplace=True)

hd_vars = pd.read_csv('/home/camille/PoolImpHuman/data/omniexpress-isec-LDHD/HDonly.coords',
                         sep='\t',
                         header=None,
                         names=['CHROM', 'POS'])
hd_vars['variants'] = pd.Series(['20:{}'.format(pos) for pos in hd_vars['POS']], dtype=str)
hd_vars.set_index('variants', inplace=True)


# Analysis of the genetic map used for linkage

pd.set_option('precision', 15)
gmap = pd.read_csv('/home/camille/phaser-factory/examples/5_snps_interpolated_HapMap2_map_20',
                   sep=' ',
                   header=None,
                   names=['ID', 'POS', 'dist'])
# modifies marker ID to chrom:pos identifier
gmap['variants'] = pd.Series(['20:{}'.format(pos) for pos in gmap['POS']], dtype=str)
gmap.set_index('variants', inplace=True)


# Lists of HD/LD markers (56,627 and 18,131 variants) (bead-chips not intersected with 1KGP data

ld_illumina = pd.read_csv('/home/camille/PoolImpHuman/data/omniexpress24/InfiniumOmniExpress-chr20-CHROM-POS.txt',
                          sep='\t',
                          header=None,
                          names=['CHROM', 'POS'])
ld_illumina['variants'] = [':'.join([chrom, pos]) for chrom, pos in zip(ld_illumina['CHROM'].astype(str).values, ld_illumina['POS'].astype(str).values)]
ld_illumina.set_index('variants', inplace=True)

hd_illumina = pd.read_csv('/home/camille/PoolImpHuman/data/omniexpress8/InfiniumOmniExpress-chr20-CHROM-POS.txt',
                          sep='\t',
                          header=None,
                          names=['CHROM', 'POS'])
hd_illumina['variants'] = [':'.join([chrom, pos]) for chrom, pos in zip(hd_illumina['CHROM'].astype(str).values, hd_illumina['POS'].astype(str).values)]
hd_illumina.set_index('variants', inplace=True)



# Discordance at markers

diff1 = quality1gt.diff().dropna()
diff1.columns = ['discordance']
if imputed_gtgp2 is not None:
     diff2 = quality2gt.diff().dropna()
     diff2.columns = ['discordance']

# Correctly imputed markers
diff1_ok = diff1.query('discordance == 0')
assert diff1_ok.mean().values == 0.0

# Correctly decoded/assayed markers
known_markers1 = dfpooled.missing_rate.query('missing_rate == 0.0')
assert known_markers1.mean().values == 0.0

# Unassayed markers
unknown_markers1 = dfpooled.missing_rate.query('missing_rate == 1.0')

# Correctly imputed and correctly decoded markers or fully assayed markers (non missing markers before imputation)
known_markers_diff1_ok = diff1.discordance.to_frame().join(dfpooled.missing_rate, how='inner').query('missing_rate == 0.0 and discordance == 0')

# Correctly imputed and unassayed markers (missing markers before imputation)
unknown_markers_diff1_ok = diff1.discordance.to_frame().join(dfpooled.missing_rate, how='inner').query('missing_rate == 1.0 and discordance == 0')

# Wrongly imputed and correctly decoded markers or fully assayed markers (non missing markers before imputation)
known_markers_diff1_ko = diff1.discordance.to_frame().join(dfpooled.missing_rate, how='inner').query('missing_rate == 0.0 and discordance == 1')

# Wrongly imputed and unassayed markers (missing markers before imputation)
unknown_markers_diff1_ko = diff1.discordance.to_frame().join(dfpooled.missing_rate, how='inner').query('missing_rate == 1.0 and discordance == 1')

# LD coordinates x Correctly decoded markers or fully assayed LD markers (non missing markers before imputation)
markers_diff1_ld = diff1.discordance.to_frame().join(dfpooled.missing_rate, how='inner').loc[ld_vars.index]
assayed_impok_ld = markers_diff1_ld.query('missing_rate == 0.0 and discordance == 0')
unassayed_impok_ld = markers_diff1_ld.query('missing_rate == 1.0 and discordance == 0')
assayed_impko_ld = markers_diff1_ld.query('missing_rate == 0.0 and discordance != 0')
unassayed_impko_ld = markers_diff1_ld.query('missing_rate == 1.0 and discordance != 0')
try:
    known_markers_diff1_ld = diff1.discordance.to_frame().join(known_markers1, how='inner').loc[ld_vars.index]
except KeyError:
    known_markers_diff1_ld = pd.DataFrame([])

# HD coordinates x Correctly decoded markers or fully assayed LD markers (non missing markers before imputation)
markers_diff1_hd = diff1.discordance.to_frame().join(dfpooled.missing_rate, how='inner').loc[hd_vars.index]
assayed_impok_hd = markers_diff1_hd.query('missing_rate == 0.0 and discordance == 0')
unassayed_impok_hd = markers_diff1_hd.query('missing_rate == 1.0 and discordance == 0')
assayed_impko_hd = markers_diff1_hd.query('missing_rate == 0.0 and discordance != 0')
unassayed_impko_hd = markers_diff1_hd.query('missing_rate == 1.0 and discordance != 0')
try:
    known_markers_diff1_hd = diff1.discordance.to_frame().join(known_markers1, how='inner').loc[hd_vars.index]
except KeyError:
    known_markers_diff1_hd = pd.DataFrame([])



dict_counts = {'true':
                   {'LD': ld_vars.shape[0], 'HD': hd_vars.shape[0]},
               'partially_assayed':
                   {'assayed':  # known_markers1.shape[0]
                        {'LD':  {'imputed': {'correct': assayed_impok_ld.shape[0],
                                             'wrong': assayed_impko_ld.shape[0]}},
                         'HD': {'imputed': {'correct': assayed_impok_hd.shape[0],
                                            'wrong': assayed_impko_hd.shape[0]}}},
                    'unassayed':  # unknown_markers1.shape[0]
                        {'LD': {'imputed': {'correct': unassayed_impok_ld.shape[0],
                                            'wrong': unassayed_impko_ld.shape[0]}},
                         'HD': {'imputed': {'correct': unassayed_impok_hd.shape[0],
                                            'wrong': unassayed_impko_hd.shape[0]}}}
                    }
               }

counts_ld = [
    ['LD', 1, 1, assayed_impok_ld.shape[0]],
    ['LD', 1, 0, assayed_impko_ld.shape[0]],
    ['LD', 0, 0, unassayed_impko_ld.shape[0]],
    ['LD', 0, 1, unassayed_impok_ld.shape[0]],
]
df_counts_ld = pd.DataFrame(counts_ld, columns=['density', 'assayed', 'imputed', 'counts'])


counts_hd = [
    ['HD', 1, 1, assayed_impok_hd.shape[0]],
    ['HD', 1, 0, assayed_impko_hd.shape[0]],
    ['HD', 0, 0, unassayed_impko_hd.shape[0]],
    ['HD', 0, 1, unassayed_impok_hd.shape[0]],
]
df_counts_hd = pd.DataFrame(counts_hd, columns=['density', 'assayed', 'imputed', 'counts'])

df_counts = pd.concat([df_counts_ld, df_counts_hd], axis=0, ignore_index=True)
df_counts.to_latex(buf=os.path.join(outdir, 'assayed_imputed_counts.tex'),
                   sparsify=True,
                   index=False,
                   caption=os.path.basename(outdir))

if False:
    # lists of HD/LD markers

    dfhdgt = vcfdf.PandasMixedVCF(truegt, format='GT', indextype='chrom:pos')
    dfpooledhd = vcfdf.PandasMixedVCF(pooledgl, format='GL', indextype='chrom:pos')

    decoded_ok = dfpooledhd.missing_rate[dfpooledhd.missing_rate == 0].dropna()
    missing_pooled = dfpooledhd.missing_rate
    af_info_pooled = dfpooledhd.af_info

    dfhd = gmap.join([af_info_pooled, missing_pooled], how='inner')
    mydiff = dfhd['dist'].diff(periods=1)

    # Phaser commit c4be822e
    # dfhd['linkage'] = mydiff
    # dfhd['linkage'].iloc[0] = 0.01  # as in VcfUtils.cpp l.512

    # Phaser commit 32c7b59e
    dfhd['linkage'] = mydiff.dropna().append(pd.Series(0.0))

    dfhd['imp-modif'] = dfhd.join(known_markers_diff2, how='left')['HG01063']
    dfhd['imp-modif'].replace(0.0, np.nan, inplace=True)

    # log the small linkage values?
    mylog10 = np.vectorize(lambda x: np.log10(x) if x > 1e-10 else -10.0)
    dfhd['log-linkage'] = mylog10(mydiff.values)

    # Where on the genetic map are the decoded markers modified by imputation?
    modified_markers = dfhd[~dfhd['imp-modif'].isna()]
    modified_markers['log-linkage'] = dfhd['log-linkage'].loc[modified_markers.index]

    ld_modified_markers = modified_markers[['ID', 'imp-modif', 'log-linkage']].join(ld_illumina, how='inner', rsuffix='_omnimap')
    hd_modified_markers = modified_markers[['ID', 'imp-modif', 'log-linkage']].join(hd_illumina, how='inner', rsuffix='_omnimap')


    # Naive simple plots

    # what is the relationship between AF_INFO in PAN reference panel and linkage between the target markers on the chosen genetic map?
    dfhd.plot.scatter(x='af_info', y='linkage', s=8)
    plt.show()

    plt.clf()
    sns.scatterplot(data=dfhd, x='af_info', y='log-linkage', hue='missing_rate', s=8)
    # add to the plot the 24 markers correctly decoded that were modified then during imputation?
    sns.scatterplot(data=dfhd[~dfhd['imp-modif'].isna()], x='af_info', y='log-linkage', color='r', s=30, legend=False)
    plt.show()


    # Plot linkage with respect to the positions of the modified markers in the genome?
    # with recombination hotspots where very high linkage too?

    plt.clf()
    sns.scatterplot(data=dfhd, x='POS', y='log-linkage', color='w', alpha=0.0, s=8)
    # high linkage regions only
    sns.scatterplot(data=dfhd[dfhd['log-linkage'] < -3.0], x='POS', y='log-linkage', color='g', s=15)
    # modified markers only
    sns.scatterplot(data=ld_modified_markers, x='POS', y='log-linkage', color='r', marker='o', s=50, label='ld_modified_markers', legend=False)
    sns.scatterplot(data=hd_modified_markers, x='POS', y='log-linkage', color='r', marker='v', s=50, label='hd_modified_markers', legend=False)
    for marker in dfhd[~dfhd['imp-modif'].isna()][['POS', 'log-linkage']].itertuples():
        plt.axvline(marker[1], color='r', lw=0.5)
    plt.legend()
    plt.show()



"""
HG01063	PUR			1
Has Affy 6.0 Genotypes

PUR	Puerto Rican	Puerto Rican in Puerto Rico	18.4	-66.1	AMR	American Ancestry	#710027	2	150	1000 Genomes on GRCh38,1000 Genomes 30x on GRCh38,The Human Genome Structural Variation Consortium,1000 Genomes phase 3 release,1000 Genomes phase 1 release

population specific af for these markers? 

20	23266616	rs1541230	T	C	100	PASS	AC=36;  AF=0.0690895;   AN=480;NS=2504;DP=17780;EAS_AF=0.0645;   AMR_AF=0.0202;  AFR_AF=0.1891;EUR_AF=0.001;SAS_AF=0.0164;AA=T|||;VT=SNP
20	54059419	rs212595	A	G	100	PASS	AC=81;  AF=0.180312;    AN=480;NS=2504;DP=13695;EAS_AF=0.1835;    AMR_AF=0.2997;  AFR_AF=0.0961;EUR_AF=0.2276;SAS_AF=0.1575;AA=A|||;VT=SNP
20	57109200	rs11699914	G	A	100	PASS	AC=17;  AF=0.0351438;   AN=480;NS=2504;DP=19853;EAS_AF=0.002;   AMR_AF=0.0519;  AFR_AF=0.0076;EUR_AF=0.0795;SAS_AF=0.0491;AA=G|||;VT=SNP
20	61714591	rs6010838	T	C	100	PASS	AC=45;  AF=0.11262; AN=480;NS=2504;DP=18667;EAS_AF=0.005;   AMR_AF=0.0893;  AFR_AF=0.261;EUR_AF=0.1123;SAS_AF=0.0399;AA=C|||;VT=SNP
"""

"""
ld_modified_markers_true_gl
Out[22]: 
                     HG01063
variants                    
20:23266616  (0.0, 0.0, 1.0)
20:54059419  (0.0, 0.0, 1.0)
20:57109200  (0.0, 0.0, 1.0)
20:61714591  (0.0, 0.0, 1.0)

ld_modified_markers_phaser_gl
Out[24]: 
                     HG01063
variants                    
20:23266616  (0.0, 0.0, 1.0)
20:54059419  (0.0, 0.0, 1.0)
20:57109200  (0.0, 0.0, 1.0)
20:61714591  (0.0, 0.0, 1.0)

ld_modified_markers_beagle_gl
Out[25]: 
                     HG01063
variants                    
20:23266616  (0.0, 0.0, 1.0)
20:54059419  (0.0, 0.0, 1.0)
20:57109200  (0.0, 0.0, 1.0)
20:61714591  (0.0, 0.0, 1.0)

# from pooled HD data
qphaser2gl.imputedobj.genotypes().loc[ld_modified_markers.index]
Out[3]: 
                                                      HG01063
variants                                                     
20:23266616   (0.0, 0.9513159990310669, 0.048684000968933105)
20:54059419  (0.0, 0.9988980293273926, 0.0011020000092685223)
20:57109200  (0.0, 0.9988809823989868, 0.0011190000222995877)
20:61714591   (0.0, 0.9913960099220276, 0.008604000322520733)


dfpooledhd.genotypes().loc[ld_modified_markers.index]
Out[5]: 
                                                       HG01063
variants                                                      
20:23266616  (-12.0, -0.15076099336147308, -0.5326970219612...
20:54059419  (-12.0, -0.2254219949245453, -0.3926349878311157)
20:57109200  (-12.0, -0.15076099336147308, -0.5326970219612...
20:61714591  (-12.0, -0.15076099336147308, -0.5326970219612...

# So, actually these markers are half-decoded only!!!
Ex.
np.power(10, (-12.0, -0.2254219949245453, -0.3926349878311157))
Out[6]: array([1.00000000e-12, 5.95083632e-01, 4.04916069e-01])

can correspond to pooling patterns:
* 3  1	0	2	2	0	1	1
* 2	 2	0	3	1	0	1	1
Looks sensible


After fixing the trinary encoding
trinary_encoding().loc[ld_modified_markers.index]
Out[47]: 
             HG01063
variants            
20:23266616       -1
20:54059419       -1
20:57109200       -1
20:61714591       -1
"""


"""
1 - dfpooledhd.missing_rate.mean()
Out[2]: 
missing_rate    0.53697553940452

len(dfpooledhd.missing_rate) - dfpooledhd.missing_rate.sum()
Out[3]: 
missing_rate    28297.0
dtype: float64

dfpooledhd.missing_rate[dfpooledhd.missing_rate == 0].dropna()
Out[5]: 
             missing_rate
variants                 
20:61651              0.0
20:63231              0.0
20:69408              0.0
20:71093              0.0
20:80655              0.0
...                   ...
20:62905851           0.0
20:62906157           0.0
20:62910548           0.0
20:62915126           0.0
20:62947458           0.0
[28297 rows x 1 columns]

"""