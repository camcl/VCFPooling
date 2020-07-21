import matplotlib.pyplot as plt

from scripts.VCFPooling.poolSNPs import dataframe as vcfdf
from scripts.VCFPooling.poolSNPs import parameters as prm
from persotools.files import *

"""
Plot AAF values after imputation (computed with Beagle)
vs. AAF values in the original dataset (theoretical).
Compare effects of imputation:
- GT serves as baseline for comparison
- adaptive GL values for filling in missing genotypes
"""

print('Configure working directory'.ljust(80, '.'))
# chunk10000_20190725 settings gave the best results for Beagle
dirs = {'default': os.path.join(prm.WD, 'gt', 'stratified', 'all_snps_all_samples'),
        'adaptive': os.path.join(prm.WD, 'gl', 'gl_adaptive', 'chunk10000_20190725')}
os.chdir(dirs['adaptive'])

print('Load files path locations'.ljust(80, '.'))
paths = {'preimp_gt_default': os.path.join(dirs['default'], prm.POOLED['b1']) + '.vcf.gz',
         'postimp_gt_default': os.path.join(dirs['default'], prm.POOLED['gtonly']) + '.vcf.gz',
         'preimp_gl_adaptive': os.path.join(dirs['adaptive'], prm.POOLED['b1']) + '.vcf.gz',
         'postimp_gl_adaptive': os.path.join(dirs['adaptive'], prm.POOLED['gtonly']) + '.vcf.gz'
         }

print('Concatenate together AAF from files to compare'.ljust(80, '.'))
basefile = os.path.join(prm.PATH_GT_FILES, 'IMP.chr20.pooled.snps.gt.chunk{}.vcf.gz'.format(prm.CHK_SZ))
pdvcfbase = vcfdf.PandasVCF(basefile, indextype='id')
afinfo = pdvcfbase.af_info.to_frame()

compaafs = vcfdf.PanelVCF(**paths).join(idt='id')

dfplot = afinfo.join(compaafs, how='inner')
dfplot.sort_values(by='af_info', inplace=True)
print('dfplot\n', dfplot)

print('Plotting'.ljust(80, '.'))
plt.rcParams["figure.figsize"] = [12, 6]
plt.rcParams["figure.autolayout"] = True
colors = ['tab:green', 'tab:olive', 'tab:cyan', 'tab:blue']

fig, axis = plt.subplots()
for k, y_key in enumerate(paths.keys()):
    dfplot.plot.scatter(x='af_info',
                        y='aaf_' + y_key,
                        ax=axis,
                        label=y_key,
                        marker='o',
                        color=colors[k],
                        s=0.7)
plt.xlabel('Theoretical AAF')
plt.ylabel('Dataset AAF')
plt.plot(range(2), linestyle='-', color='k')
plt.legend()
plt.title('AAF evolution through processing of data sets', loc='center')
plt.suptitle("")
# plt.savefig('aaf_evol.scatter.gtgl.pooled.chunk{}.png'.format(prm.CHK_SZ),
#             orientation='landscape',
#             dpi='figure')
plt.show()

###
# df_err = pd.read_csv('pooled.sorted.chunk{}.csv',
#                      sep='\t',
#                      encoding='utf-8',
#                      index_col=0,
#                      usecols=[0] + list(range(3, prm.NB_IMP + 3)),
#                      skiprows=[0, 1]
#                      )
# allplt.plot_aaf_twist('pooled',
#                       df_err,
#                       os.path.join(dirs['adaptive'], 'IMP.chr20.pooled.imputed.gt.chunk{}.vcf.gz'.format(prm.CHK_SZ)))

