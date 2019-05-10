import subprocess
import numpy as np

from scripts.poolSNPs import parameters as prm
from scripts.poolSNPs import patterns as pat
from scripts.poolSNPs import bcftools
from scripts.poolSNPs.alleles import alleles_tools as alltls
from persotools.debugging import *
from persotools.files import *

'''
Run bcftools manipulations for preprocessing the files which are used for Beagle imputation.
Run Beagle.
'''
chk_sz = prm.CHK_SZ
pardir = os.path.join(prm.WD, prm.GTGL.lower())
if prm.GTGL == 'GT':
    #cd = pardir
    cd = os.path.join(pardir, 'stratified')
if prm.GTGL == 'GL':
    if prm.unknown_gl != 'adaptative':
        cd = os.path.join(pardir, 'gl_' + '_'.join(np.around(prm.unknown_gl, decimals=2).astype(str)))
    else:
        cd = os.path.join(pardir, 'gl_adaptative')
    mkdir(cd)
os.chdir(cd)

raw = prm.RAW
pooled = prm.POOLED
missing = prm.MISSING
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
    print('{} compressed to {}'.format(os.path.join(pardir, dic['vcf']),
                                       dic['gz']))
    if os.path.exists(os.path.join(pardir, dic['vcf'])):
        bcftools.bgzip(os.path.join(pardir, dic['vcf']),
                       dic['gz'],
                       cd) # bgzip the file in the corresponding GL (GT) folder for missing values
        bcftools.sort(dic['gz'], cd)
        bcftools.index(dic['gz'], cd)
        #delete_file(pardir + dic['vcf'])
    if os.path.exists(os.path.join(cd, dic['vcf'])):
        bcftools.bgzip(os.path.join(cd, dic['vcf']),
                       dic['gz'],
                       cd) # bgzip the file in the corresponding GL (GT) folder for missing values
        bcftools.sort(dic['gz'], cd)
        bcftools.index(dic['gz'], cd)
        #delete_file(pardir + dic['vcf'])

print('Set size for REF: ', prm.NB_REF)
print('Set size for IMP: ', prm.NB_IMP)
samples_files = ['cat {}/ALL.chr20.snps.allID.txt '.format(pardir)
                 + '| head -{} > {}/ALL.chr20.snps.impID.txt'.format(prm.NB_IMP, pardir),
                 'cat {}/ALL.chr20.snps.allID.txt '.format(pardir)
                 + '| tail -{} > {}/ALL.chr20.snps.refID.txt'.format(prm.NB_REF, pardir),
                 'dos2unix {}/ALL.chr20.snps.refID.txt'.format(pardir),
                 'dos2unix {}/ALL.chr20.snps.impID.txt'.format(pardir)]
for f in samples_files:
    subprocess.run(f, shell=True, cwd=cd)

### REF/IMP SAMPLING
print('\n\nREF/IMP SAMPLING'.ljust(80, '.'))
for (folder, dic) in tuple([(cd, pooled),
                            (cd, missing),
                            (prm.WD + '/gt/stratified/', raw)]): # cd if gt stratified
    delete_file(os.path.join(folder, dic['imp']))
    delete_file(os.path.join(folder, dic['imp'] + '.csi'))
    bcftools.sampling(dic['gz'],
                      dic['imp'],
                      '{}/ALL.chr20.snps.impID.txt'.format(pardir),
                      folder)
    bcftools.index(dic['imp'], folder)

bcftools.sampling(raw['gz'],
                  raw['ref'],
                  '{}/ALL.chr20.snps.refID.txt'.format(pardir),
                  prm.WD + '/gt/stratified')
bcftools.index(raw['ref'], cd) # prm.WD + '/gt/'

for (folder, f) in tuple([(prm.WD + '/gt/stratified', raw['imp']),
                          (prm.WD + '/gt/stratified', raw['ref']),
                          (cd, pooled['imp']),
                          (cd, missing['imp'])]):
    bcftools.sort(f, folder)
    bcftools.index(f, folder)

if GTGL == 'GT':
    # TODO: test that piece of code
    folder = prm.WD + '/gt/stratified'
    bcftools.sampling(pooled['gz'],
                      pooled['ref'],
                      '{}/ALL.chr20.snps.refID.txt'.format(pardir),
                      folder)
    bcftools.index(pooled['ref'], folder)

### GL CONVERSION FOR IMP/REF: should be done from the GT files, not possible by splitting the converted ALL.gl file
if GTGL == 'GL':
    print('GT to GL in {}'.format(os.getcwd()).ljust(80, '.'))
    alltls.file_likelihood_converter(os.path.join(prm.WD,
                                                  'gt',
                                                  raw['ref'].replace('.gl', '.gt')),
                                     raw['ref'][:-3])
    delete_file(raw['ref'])
    delete_file(raw['ref'] + '.csi')
    bcftools.bgzip(raw['ref'][:-3], raw['ref'], cd)
    bcftools.sort(raw['ref'], cd)
    bcftools.index(raw['ref'], cd)
    delete_file(raw['ref'][:-3])

    # for dic in [pooled, missing]: # raw
    #     if prm.unknown_gl == 'adaptative':
    #         pat.adaptative_likelihood_converter(os.path.join(prm.WD,
    #                                                          'gt',
    #                                                          dic['imp'].replace('.gl', '.gt')),
    #                                             dic['imp'][:-3])
    #     else:
    #         alltls.file_likelihood_converter(os.path.join(prm.WD,
    #                                                       'gt',
    #                                                       dic['imp'].replace('.gl', '.gt')),
    #                                          dic['imp'][:-3])
    #
    #     delete_file(dic['imp'])
    #     delete_file(dic['imp'] + '.csi')
    #     bcftools.bgzip(dic['imp'][:-3], dic['imp'], cd)
    #     bcftools.sort(dic['imp'], cd)
    #     bcftools.index(dic['imp'], cd)
    #     delete_file(dic['imp'][:-3])

### BEAGLE ROUND#1: PHASING
print('\n\nBEAGLE ROUND#1'.ljust(80, '.'))
# REF dataset is GT formatted (if directly available)
delete_file(raw['b1r'] + '.vcf.gz')
delete_file(raw['b1i'] + '.vcf.gz')

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
                   + os.path.join(prm.WD, 'gt/stratified', raw['imp'].replace('.gl', '.gt')),
                   'impute=false',
                   'gprobs=true',
                   'out=' + raw['b1i'],
                   '&',
                   'java -Xmx4000m -jar {}'.format(prm.BEAGLE_JAR),
                   '{}='.format('gt')
                   + os.path.join(prm.WD, 'gt/stratified', raw['ref'].replace('.gl', '.gt')),
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
bcftools.index(raw['b1r'] + '.vcf.gz', cd)
bcftools.index(raw['b1i'] + '.vcf.gz', cd)
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
    bcftools.index(dic['b1'] + '.vcf.gz', cd)
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
    bcftools.index(dic['cfgt'] + '.vcf.gz', cd)


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
    bcftools.index(dic['b2'] + '.vcf.gz', cd)


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
    bcftools.index(dic['corr'] + '.vcf.gz', cd)

    gtonly = ' '.join(["bcftools annotate -x 'FORMAT'",
                       dic['corr'] + '.vcf.gz',
                       "| bgzip >",
                       dic['gtonly'] + '.vcf.gz'])

    subprocess.run(gtonly, shell=True, cwd=cd)
    bcftools.index(dic['gtonly'] + '.vcf.gz', cd)


### REMOVE .log FILES and cfgt created
for f in os.scandir(cd):
    if f.is_file()and (f.path.endswith('.log') or '.cfgt.' in f.path):
        delete_file(f)
