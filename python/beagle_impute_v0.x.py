import subprocess
import numpy as np

from scripts.VCFPooling.poolSNPs import parameters as prm
from scripts.VCFPooling.poolSNPs import pybcf
from scripts.VCFPooling.python.archived.alleles import alleles_tools as alltls

from persotools.files import *
from persotools.struct import *

'''
Run bcftools manipulations for preprocessing the files which are used for Beagle imputation.
Run Beagle.
'''

### PARAMETERS
chk_sz = prm.CHK_SZ
path_gt_files = prm.PATH_GT_FILES
path_gl_files = prm.PATH_GL_FILES

if prm.GTGL == 'GT':
    cd = path_gt_files
if prm.GTGL == 'GL':
    if prm.unknown_gl != 'adaptive':
        cd = os.path.join(path_gl_files, 'gl_' + '_'.join(np.around(prm.unknown_gl, decimals=2).astype(str)))
    else:
        cd = os.path.join(path_gl_files, 'gl_adaptive')
    mkdir(cd)
os.chdir(cd)

raw = NamedDict('raw', list(prm.RAW.keys()), list(prm.RAW.values()))
pooled = NamedDict('pooled', list(prm.POOLED.keys()), list(prm.POOLED.values()))
missing = NamedDict('missing', list(prm.MISSING.keys()), list(prm.MISSING.values()))
subset = prm.SUBSET #if susbset get the first 1000 lines to recreate ech method file: pooled, missing etc
trc = prm.SUBCHUNK
GTGL = prm.GTGL

print('Start processing files for running BEAGLE'.ljust(80, '.'))

if subset:
    print('Subset main set in {}'.format(os.getcwd()).ljust(80, '.'))
    for dic in [pooled, missing, raw]:
        for k, v in dic.items():
            dic[k] = v.replace('chunk' + str(chk_sz), 'chunk' + str(trc))

### BGZIP ALL
print('\n\nBGZIP in {}'.format(os.getcwd()).ljust(80, '.'))
for dic in [pooled, missing]:

    if prm.GTGL == 'GT':
        delete_file(dic['gz'] + '.csi')
        print('{} compressed to {}'.format(os.path.join(path_gt_files, dic['vcf']),
                                           os.path.join(cd, dic['gz'])))
        pybcf.bgzip(os.path.join(path_gt_files, dic['vcf']),
                    dic['gz'],
                    cd) # bgzip the file in the corresponding GT folder for missing values
        pybcf.sort(dic['gz'], cd)
        pybcf.index(dic['gz'], cd)
        #delete_file(path_gt_files + dic['vcf'])

    if prm.GTGL == 'GL':
        delete_file(dic['gz'] + '.csi')
        print('{} compressed to {}'.format(os.path.join(path_gl_files, dic['vcf']),
                                           os.path.join(cd, dic['gz'])))
        pybcf.bgzip(os.path.join(path_gl_files, dic['vcf']),
                    dic['gz'],
                    cd) # bgzip the file in the corresponding GL folder for missing values
        pybcf.sort(dic['gz'], cd)
        pybcf.index(dic['gz'], cd)
        #delete_file(path_gl_files + dic['vcf'])

print('Set size for REF (reference panel): ', prm.NB_REF)
print('Set size for IMP (study population): ', prm.NB_IMP)
samples_files = ['cat {}/ALL.chr20.snps.allID.txt '.format(prm.WD + '/gt')
                 + '| head -{} > {}/ALL.chr20.snps.impID.txt'.format(prm.NB_IMP, prm.WD + '/gt'),
                 'cat {}/ALL.chr20.snps.allID.txt '.format(prm.WD + '/gt')
                 + '| tail -{} > {}/ALL.chr20.snps.refID.txt'.format(sum([prm.NB_REF, prm.NB_IMP]), prm.WD + '/gt'),
                 'dos2unix {}/ALL.chr20.snps.refID.txt'.format(prm.WD + '/gt'),
                 'dos2unix {}/ALL.chr20.snps.impID.txt'.format(prm.WD + '/gt')]
for f in samples_files:
    subprocess.run(f, shell=True, cwd=cd)

### REF/IMP SAMPLING
print('\n\nREF/IMP SAMPLING'.ljust(80, '.'))
for (folder, dic) in tuple([(cd, pooled),
                            (cd, missing),
                            (path_gt_files, raw)]):
    delete_file(os.path.join(folder, dic['imp']))
    delete_file(os.path.join(folder, dic['imp'] + '.csi'))
    pybcf.sampling(dic['gz'],
                   dic['imp'],
                      '{}/ALL.chr20.snps.impID.txt'.format(prm.WD + '/gt'),
                   folder)
    if dic != raw: #TODO: move aftewr beagle round1
        pybcf.rename_samples(dic['imp'], dic['imp'][:-3], cd, '_IMP')
        pybcf.bgzip(dic['imp'][:-3], dic['imp'], cd)
        delete_file(dic['imp'][:-3])
    # remove str(var.POS) == '59973567'
    if False:
        cmd = ' '.join(['bcftools view -e POS=59973567',
                        '-Oz -o',
                        dic['imp'],
                        dic['imp']])
        subprocess.run(cmd, shell=True, cwd=folder)
    pybcf.index(dic['imp'], folder)

pybcf.sampling(raw['gz'],
               raw['ref'],
                  '{}/ALL.chr20.snps.refID.txt'.format(prm.WD + '/gt'),
               path_gt_files)
pybcf.index(raw['ref'], path_gt_files)

for (folder, f) in tuple([(path_gt_files, raw['imp']),
                          (path_gt_files, raw['ref']),
                          (cd, pooled['imp']),
                          (cd, missing['imp'])]):
    pybcf.sort(f, folder)
    pybcf.index(f, folder)

if GTGL == 'GT':
    folder = path_gt_files
    pybcf.sampling(pooled['gz'],
                   pooled['ref'],
                      '{}/ALL.chr20.snps.refID.txt'.format(prm.WD + '/gt'),
                   folder)
    pybcf.index(pooled['ref'], folder)

### GL CONVERSION FOR IMP/REF: should be done from the GT files, not possible by splitting the converted ALL.gl file
if GTGL == 'GL':
    print('GT to GL in {}'.format(os.getcwd()).ljust(80, '.'))
    alltls.file_likelihood_converter(os.path.join(path_gt_files,
                                                  raw['ref'].replace('.gl', '.gt')),
                                     raw['ref'][:-3])
    delete_file(raw['ref'])
    delete_file(raw['ref'] + '.csi')
    pybcf.bgzip(raw['ref'][:-3], raw['ref'], cd)
    pybcf.sort(raw['ref'], cd)
    pybcf.index(raw['ref'], cd)
    delete_file(raw['ref'][:-3])

### BEAGLE ROUND#1: PHASING
print('\n\nBEAGLE ROUND#1'.ljust(80, '.'))
# if reference panel is GT formatted
delete_file(raw['b1r'] + '.vcf.gz')
delete_file(raw['b1i'] + '.vcf.gz')

# elif reference panel is GL encoded:
# bgl1gtgl = ' '.join(['java -Xmx4000m -jar {}'.format(prm.BEAGLE_JAR),
#                      '{}='.format('gtgl') + raw['imp'],
#                      'impute=false',
#                      'gprobs=true',
#                      'out=' + 'temp.' + raw['imp'][:-7],
#                      '&',
#                      'java -Xmx4000m -jar {}'.format(prm.BEAGLE_JAR),
#                      '{}='.format('gtgl') + raw['ref'],
#                      'impute=false',
#                      'gprobs=true',
#                      'out=' + 'temp.' + raw['ref'][:-7]
#                      ])
#
# bgl1gt = ' '.join(['java -Xmx4000m -jar {}'.format(prm.BEAGLE_JAR),
#                    '{}='.format('gt')
#                    + '{}'.format('temp.' if GTGL == 'GL' else '')
#                    + raw['imp'],
#                    'impute=false',
#                    'gprobs=true',
#                    'out=' + raw['b1i'],
#                    '&',
#                    'java -Xmx4000m -jar {}'.format(prm.BEAGLE_JAR),
#                    '{}='.format('gt')
#                    + '{}'.format('temp.' if GTGL == 'GL' else '')
#                    + raw['ref'],
#                    'impute=false',
#                    'gprobs=true',
#                    'out=' + raw['b1r']
#                    ])

bgl1gt = ' '.join(['java -Xmx4000m -jar {}'.format(prm.BEAGLE_JAR),
                   '{}='.format('gt')
                   + os.path.join(path_gt_files, raw['imp'].replace('.gl', '.gt')),
                   'impute=false',
                   'gprobs=true',
                   'out=' + raw['b1i'],
                   '&',
                   'java -Xmx4000m -jar {}'.format(prm.BEAGLE_JAR),
                   '{}='.format('gt')
                   + os.path.join(path_gt_files, raw['ref'].replace('.gl', '.gt')),
                   'impute=false',
                   'gprobs=true',
                   'out=' + raw['b1r']
                   ])

idxb1 = ' '.join(['bcftools',
                  'index -f',
                  raw['b1r'] + '.vcf.gz',
                  '&',
                  'bcftools',
                  'index -f',
                  raw['b1i'] + '.vcf.gz'
                  ])
# if GTGL == 'GL':
#     subprocess.run(bgl1gtgl, shell=True, cwd=cd)
subprocess.run(bgl1gt, shell=True, cwd=cd)
pybcf.index(raw['b1r'] + '.vcf.gz', cd)
pybcf.index(raw['b1i'] + '.vcf.gz', cd)
delete_file('temp.' + raw['imp'])
delete_file('temp.' + raw['ref'])

for dic in [pooled, missing]:
    delete_file(dic['b1'] + '.vcf.gz')

    bgl1gtgl = ' '.join(['java -Xmx4000m -jar {}'.format(prm.BEAGLE_JAR),
                         '{}='.format('gtgl') + dic['imp'],
                         'impute=false',
                         'gprobs=true',
                         'out=' + 'temp.b1'
                         ])

    bgl1gt = ' '.join(['java -Xmx4000m -jar {}'.format(prm.BEAGLE_JAR),
                       '{}='.format('gt')
                       + '{}'.format('temp.b1.vcf.gz' if GTGL == 'GL' else dic['imp']),
                       'impute=false',
                       'gprobs=true',
                       'out=' + dic['b1']
                       ])

    if GTGL == 'GL':
        subprocess.run(bgl1gtgl, shell=True, cwd=cd)
    subprocess.run(bgl1gt, shell=True, cwd=cd)
    pybcf.index(dic['b1'] + '.vcf.gz', cd)
    delete_file('temp.b1' + '.vcf.gz')

### CONFORM-GT
# conform-gt algo should be used if REF and IMP datasets contain different sets of markers
# chrom-POS: min and max from the original file
print('\n\nCONFORM-GT'.ljust(80, '.'))
# GT for reference files, even when working with GLs
for dic in [pooled, missing]:
    delete_file(dic['cfgt'] + '.vcf.gz')

    cfgt = ' '.join(['java -jar {}'.format(prm.CFGT_JAR),
                     '{}='.format('gt') + dic['b1'] + '.vcf.gz',
                     'chrom=20:60343-62965354',
                     'ref={}'.format(raw['b1r'] + '.vcf.gz'),
                     'out=' + dic['cfgt']
                     ])

    subprocess.run(cfgt, shell=True, cwd=cd)
    pybcf.index(dic['cfgt'] + '.vcf.gz', cd)


### BEAGLE (ROUND#2): IMPUTING
print('\n\nBEAGLE (ROUND#2)'.ljust(80, '.'))
for dic in [pooled, missing]:
    delete_file(dic['b2'] + '.vcf.gz')

    bgl2 = ' '.join(['java -Xmx4000m -jar {}'.format(prm.BEAGLE_JAR),
                     'gt=' + dic['cfgt'] + '.vcf.gz',
                     'ref={}.vcf.gz'.format(raw['b1r']),
                     'impute=true',
                     'gprobs=true',
                     'out=' + dic['b2']
                     ])

    subprocess.run(bgl2, shell=True, cwd=cd)
    pybcf.index(dic['b2'] + '.vcf.gz', cd)


### FIX DS AND GP FORMAT FIELDS, SORT OUT GT
print('\n\nREFORMATTING GP AND DS FIELDS'.ljust(80, '.'))
for dic in [pooled, missing]:
    delete_file(dic['corr'] + '.vcf.gz')

    refmt = ' '.join(["bcftools view {}.vcf.gz".format(dic['b2']),
                      "| sed 's/##FORMAT=<ID=DS,Number=A,Type=Float/##FORMAT=<ID=DS,Number=1,Type=String/'",
                      "| sed 's/##FORMAT=<ID=GP,Number=G,Type=Float/##FORMAT=<ID=GP,Number=3,Type=String/'",
                      "| bcftools view -Oz -o {}.vcf.gz".format(dic['corr'])
                      ])

    subprocess.run(refmt, shell=True, cwd=cd)
    pybcf.index(dic['corr'] + '.vcf.gz', cd)

    gtonly = ' '.join(["bcftools annotate -x 'FORMAT'",
                       dic['corr'] + '.vcf.gz',
                       "| bgzip >",
                       dic['gtonly'] + '.vcf.gz'])

    subprocess.run(gtonly, shell=True, cwd=cd)
    pybcf.index(dic['gtonly'] + '.vcf.gz', cd)


### REMOVE .log FILES and cfgt created
for f in os.scandir(cd):
    if f.is_file()and (f.path.endswith('.log') or '.cfgt.' in f.path):
        delete_file(f)
