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
Is there a correlation between the markers linkage and what data/how much is missing?
Look first at a single sample HG01063 (Phaser works on single samples anyways)
'''


# Data parameters and plotting features

rQ = 1000
bS = 0.01
x_data = 'binned_maf'  # 'binned_maf_info'

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


# Configure data/plots paths

outdir = '/home/camille/PoolImpHuman/results/20201109'
if not os.path.exists(outdir):
    os.mkdir(outdir)

truegt = '/home/camille/PoolImpHuman/data/20200906a/sHG01063.IMP.chr20.snps.gt.vcf.gz'
truegl = '/home/camille/PoolImpHuman/data/20200906a/sHG01063.IMP.chr20.snps.gl.vcf.gz'
pooledgl = '/home/camille/PoolImpHuman/data/20200906a/sHG01063.IMP.chr20.pooled.snps.gl.vcf.gz'
imputed_phaser1 = '/home/camille/PoolImpHuman/data/20200906a/sHG01063.IMP.chr20.pooled.snps.gl.full.postgenos.vcf.gz'
imputed_phaser2 = '/home/camille/PoolImpHuman/data/20201110/sHG01063.IMP.chr20.pooled.snps.gl.full.postgenos.vcf.gz'

qphaser1gt = qual.QualityGT(truegt, imputed_phaser1, 0, idx='chrom:pos')
qphaser1gl = qual.QualityGL(truegl, imputed_phaser1, 0, idx='chrom:pos')

qphaser2gt = qual.QualityGT(truegt, imputed_phaser2, 0, idx='chrom:pos')
qphaser2gl = qual.QualityGL(truegl, imputed_phaser2, 0, idx='chrom:pos')

dfpooled = vcfdf.PandasMixedVCF(pooledgl, format='GL', indextype='chrom:pos')

# LDHD imputation
imputed_beagle_ldhd = '/home/camille/PoolImpHuman/data/20201029/sHG01063.IMP.chr20.pooled.imputed.vcf.gz'
imputed_phaser_ldhd = '/home/camille/PoolImpHuman/data/20201031/sHG01063.IMP.chr20.snps.gl.full.postgenos.vcf.gz'

qbeagle_ldhd = qual.QualityGL(truegl, imputed_beagle_ldhd, 0, idx='chrom:pos')
qphaser_ldhd = qual.QualityGL(truegl, imputed_phaser_ldhd, 0, idx='chrom:pos')

# Discordance at markers
diff1 = qphaser1gt.diff().dropna()
diff2 = qphaser2gt.diff().dropna()

# Discordance at markers where genotype was decoded i.e. not missing
known_markers_diff1 = diff1.HG01063[(dfpooled.missing_rate == 0.0).values.flatten()]
known_markers_diff2 = diff2.HG01063[(dfpooled.missing_rate == 0.0).values.flatten()]

# Discordance at markers where genotype was not decoded i.e. missing
unknown_markers_diff1 = diff1.HG01063[(dfpooled.missing_rate == 1.0).values.flatten()]
unknown_markers_diff2 = diff2.HG01063[(dfpooled.missing_rate == 1.0).values.flatten()]

# comp = known_markers_diff1.to_frame().join(known_markers_diff2, how='inner', lsuffix='_1', rsuffix='_2')
# comp[comp.HG01063_1 == 1]

# Analysis of the genetic map used for linkage
pd.set_option('precision', 15)
gmap = pd.read_csv('/home/camille/phaser-factory/examples/5_snps_interpolated_HapMap2_map_20',
                   sep=' ',
                   header=None,
                   names=['ID', 'POS', 'dist'])
# modifies marker ID to chrom:pos identifier
gmap['variants'] = pd.Series(['20:{}'.format(pos) for pos in gmap['POS']], dtype=str)
gmap.set_index('variants', inplace=True)

# lists of HD/LD markers
ld_markers = pd.read_csv('/home/camille/PoolImpHuman/data/omniexpress24/InfiniumOmniExpress-chr20-CHROM-POS.txt',
                         sep='\t',
                         header=None,
                         names=['CHROM', 'POS'])
ld_markers['variants'] = [':'.join([chrom, pos]) for chrom, pos in zip(ld_markers['CHROM'].astype(str).values, ld_markers['POS'].astype(str).values)]
ld_markers.set_index('variants', inplace=True)

hd_markers = pd.read_csv('/home/camille/PoolImpHuman/data/omniexpress8/InfiniumOmniExpress-chr20-CHROM-POS.txt',
                         sep='\t',
                         header=None,
                         names=['CHROM', 'POS'])
hd_markers['variants'] = [':'.join([chrom, pos]) for chrom, pos in zip(hd_markers['CHROM'].astype(str).values, hd_markers['POS'].astype(str).values)]
hd_markers.set_index('variants', inplace=True)

dfhdgt = vcfdf.PandasMixedVCF(truegt, format='GT', indextype='chrom:pos')
dfpooledhd = vcfdf.PandasMixedVCF(pooledgl, format='GL', indextype='chrom:pos')

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

ld_modified_markers = modified_markers[['ID', 'imp-modif', 'log-linkage']].join(ld_markers, how='inner', rsuffix='_omnimap')
hd_modified_markers = modified_markers[['ID', 'imp-modif', 'log-linkage']].join(hd_markers, how='inner', rsuffix='_omnimap')


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


# Comparison with LDHD imputations
ld_modified_markers_true_gl = qphaser_ldhd.trueobj.genotypes().loc[ld_modified_markers.index]
ld_modified_markers_phaser_gl = qphaser_ldhd.imputedobj.genotypes().loc[ld_modified_markers.index]
ld_modified_markers_beagle_gl = qbeagle_ldhd.imputedobj.genotypes().loc[ld_modified_markers.index]



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


def trinary_encoding():
    vcfobj = dfpooledhd.load()
    vars = dfpooledhd.variants
    arr = np.empty((len(vars), len(dfpooledhd.samples)), dtype=float)
    if dfpooledhd.fmt.upper() == 'GT':
        missing = np.vectorize(lambda x: np.nan if x is None else x)
        for i, var in enumerate(vcfobj):
            # missing are read as None
            gts = np.array([g[dfpooledhd.fmt] for g in var.samples.values()]).astype(float)
            tri = missing(gts).sum(axis=-1)
            arr[i, :] = np.nan_to_num(tri, nan=-1)
    elif dfpooledhd.fmt.upper() == 'GL':
        gtnan = np.array([np.nan, np.nan, np.nan])
        gt0 = np.array([0., 0., 0.])
        gt1 = np.array([0., 1., 0.])
        gt2 = np.array([0., 0., 2.])
        for i, var in enumerate(vcfobj):
            # convert GL to trinary
            missing = lambda x: gtnan if 0.0 not in x else (gt0 if x[0] == 0.0 else (gt1 if x[1] == 0.0 else gt2))
            gts = np.array([g[dfpooledhd.fmt] for g in var.samples.values()]).astype(float)
            tri = np.apply_along_axis(missing, -1, gts).sum(axis=-1)
            arr[i, :] = np.nan_to_num(tri, nan=-1)
    dftrinary = pd.DataFrame(arr, index=vars, columns=dfpooledhd.samples, dtype=int)

    return dftrinary