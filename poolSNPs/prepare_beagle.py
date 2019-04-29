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
cd = os.path.join(prm.WD, prm.GTGL.lower())
if prm.GTGL == 'GL':
    if prm.unknown_gl != 'adaptative':
        cd = os.path.join(cd, 'gl_' + '_'.join(np.around(prm.unknown_gl, decimals=2).astype(str)))
    else:
        cd = os.path.join(cd, 'gl_adaptative')
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
if GTGL == 'GT':
    print('\n\nBGZIP in {}'.format(os.getcwd()).ljust(80, '.'))
    for dic in [pooled, missing]:
        print('{} compressed to {}'.format(dic['vcf'], dic['gz']))
        delete_file(dic['gz'] + '.csi')

        bgzip = ' '.join(['bcftools',
                          'view',
                          '-Oz -o',
                          dic['gz'],
                          dic['vcf']
                          ])

        sort = ' '.join(['bcftools',
                         'sort',
                         '-Oz -o',
                         dic['gz'],
                         dic['gz']
                         ])

        idxgz = ' '.join(['bcftools',
                          'index -f',
                          dic['gz']
                          ])

        subprocess.run(bgzip, shell=True, cwd=cd)
        subprocess.run(sort, shell=True, cwd=cd)
        subprocess.run(idxgz, shell=True, cwd=cd)
        delete_file(dic['vcf'])

    print('Set size for REF: ', prm.NB_REF)
    print('Set size for IMP: ', prm.NB_IMP)
    samples_files = ['cat ALL.chr20.snps.allID.txt | head -{} > ALL.chr20.snps.impID.txt'.format(prm.NB_IMP),
                     'cat ALL.chr20.snps.allID.txt | tail -{} > ALL.chr20.snps.refID.txt'.format(prm.NB_REF),
                     'dos2unix ALL.chr20.snps.refID.txt',
                     'dos2unix ALL.chr20.snps.impID.txt']
    for f in samples_files:
        subprocess.run(f, shell=True, cwd=cd)

### REF/IMP SAMPLING
if GTGL == 'GT':
    print('\n\nREF/IMP SAMPLING'.ljust(80, '.'))
    for dic in [pooled, missing, raw]:
        delete_file(dic['imp'])
        delete_file(dic['imp'] + '.csi')

        samp = ' '.join(['bcftools',
                         'view',
                         '-Oz -o',
                         dic['imp'],
                         '-S {}/ALL.chr20.snps.impID.txt'.format(prm.WD + '/gt'),
                         dic['gz']
                         ])

        idxsp = ' '.join(['bcftools',
                          'index -f',
                          dic['imp']
                          ])

        subprocess.run(samp, shell=True, cwd=cd)
        subprocess.run(idxsp, shell=True, cwd=cd)

    samp = ' '.join(['bcftools',
                     'view',
                     '-Oz -o',
                     raw['ref'],
                     '-S {}/ALL.chr20.snps.refID.txt'.format(prm.WD + '/gt'),
                     raw['gz']
                     ])

    idxsp = ' '.join(['bcftools',
                      '-f index',
                      raw['ref']
                      ])

    delete_file(raw['ref'] + '.csi')
    subprocess.run(samp, shell=True, cwd=cd)
    subprocess.run(idxsp, shell=True, cwd=cd)

    for f in [raw['imp'], raw['ref'], pooled['imp'], missing['imp']]:
        subprocess.run('bcftools sort {} {}'.format(f, f), shell=True, cwd=cd)

### GL CONVERSION FOR IMP/REF: should be done from the GT files, not possible by splitting the converted ALL.gl file
if GTGL == 'GL':
    print('GT to GL in {}'.format(os.getcwd()).ljust(80, '.'))
    if prm.unknown_gl == 'adaptative':
        pat.adaptative_likelihood_converter(os.path.join(prm.WD,
                                                         'gt',
                                                         raw['ref'].replace('.gl', '.gt')),
                                            raw['ref'][:-3])
    else:
        alltls.file_likelihood_converter(os.path.join(prm.WD,
                                                      'gt',
                                                      raw['ref'].replace('.gl', '.gt')),
                                         raw['ref'][:-3])

    # bgzip = ' '.join(['bcftools',
    #                   'view',
    #                   '-Oz -o',
    #                   raw['ref'],
    #                   raw['ref'][:-3]
    #                   ])
    #
    # sort = ' '.join(['bcftools',
    #                  'sort',
    #                  '-Oz -o',
    #                  raw['ref'],
    #                  raw['ref']
    #                  ])
    #
    # idxgz = ' '.join(['bcftools',
    #                   'index -f',
    #                   raw['ref']
    #                   ])
    #
    # subprocess.run(bgzip, shell=True, cwd=cd)
    # subprocess.run(sort, shell=True, cwd=cd)
    # subprocess.run(idxgz, shell=True, cwd=cd)
        
    for dic in [pooled, missing]: # raw
        if prm.unknown_gl == 'adaptative':
            pat.adaptative_likelihood_converter(os.path.join(prm.WD,
                                                             'gt',
                                                             dic['imp'].replace('.gl', '.gt')),
                                                dic['imp'][:-3])
        else:
            alltls.file_likelihood_converter(os.path.join(prm.WD,
                                                          'gt',
                                                          dic['imp'].replace('.gl', '.gt')),
                                             dic['imp'][:-3])

        delete_file(dic['imp'])
        delete_file(dic['imp'] + '.csi')

        bgzip = ' '.join(['bcftools',
                          'view',
                          '-Oz -o',
                          dic['imp'],
                          dic['imp'][:-3]
                          ])

        sort = ' '.join(['bcftools',
                         'sort',
                         '-Oz -o',
                         dic['imp'],
                         dic['imp']
                         ])

        idxgz = ' '.join(['bcftools',
                          'index -f',
                          dic['imp']
                          ])

        subprocess.run(bgzip, shell=True, cwd=cd)
        subprocess.run(sort, shell=True, cwd=cd)
        subprocess.run(idxgz, shell=True, cwd=cd)

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
                   + os.path.join(prm.WD, 'gt', raw['imp'].replace('.gl', '.gt')),
                   'impute=false',
                   'gprobs=true',
                   'out=' + raw['b1i'],
                   '&',
                   'java -Xmx4000m -jar {}'.format(prm.BEAGLE_JAR),
                   '{}='.format('gt')
                   + os.path.join(prm.WD, 'gt', raw['ref'].replace('.gl', '.gt')),
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
subprocess.run(idxb1, shell=True, cwd=cd)
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

    idxb1 = ' '.join(['bcftools',
                      'index -f',
                      dic['b1'] + '.vcf.gz'
                      ])

    if GTGL == 'GL':
        subprocess.run(bgl1gtgl, shell=True, cwd=cd)
    subprocess.run(bgl1gt, shell=True, cwd=cd)
    subprocess.run(idxb1, shell=True, cwd=cd)
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

    idxcf = ' '.join(['bcftools',
                      'index -f',
                      dic['cfgt'] + '.vcf.gz'
                      ])

    print(cfgt)

    subprocess.run(cfgt, shell=True, cwd=cd)
    subprocess.run(idxcf, shell=True, cwd=cd)


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

    idxb2 = ' '.join(['bcftools',
                      'index -f',
                      dic['b2'] + '.vcf.gz'
                      ])

    subprocess.run(bgl2, shell=True, cwd=cd)
    subprocess.run(idxb2, shell=True, cwd=cd)


### FIX DS AND GP FORMAT FIELDS, SORT OUT GT
print('\n\nREFORMATTING GP AND DS FIELDS'.ljust(80, '.'))
for dic in [pooled, missing]:
    delete_file(dic['corr'] + '.vcf.gz')

    refmt = ' '.join(["bcftools view {}.vcf.gz".format(dic['b2']),
                      "| sed 's/##FORMAT=<ID=DS,Number=A,Type=Float/##FORMAT=<ID=DS,Number=1,Type=String/'",
                      "| sed 's/##FORMAT=<ID=GP,Number=G,Type=Float/##FORMAT=<ID=GP,Number=3,Type=String/'",
                      "| bcftools view -Oz -o {}.vcf.gz".format(dic['corr'])
                      ])

    idxfm = ' '.join(['bcftools',
                      'index -f',
                      dic['corr'] + '.vcf.gz'
                      ])

    subprocess.run(refmt, shell=True, cwd=cd)
    subprocess.run(idxfm, shell=True, cwd=cd)

    gtonly = ' '.join(["bcftools annotate -x 'FORMAT'",
                       dic['corr'] + '.vcf.gz',
                       "| bgzip >",
                       dic['gtonly'] + '.vcf.gz'])

    idxgt = ' '.join(['bcftools',
                      'index -f',
                      dic['gtonly'] + '.vcf.gz'
                      ])

    subprocess.run(gtonly, shell=True, cwd=cd)
    subprocess.run(idxgt, shell=True, cwd=cd)


### REMOVE .log FILES and cfgt created
for f in os.scandir(cd):
    if f.is_file()and (f.path.endswith('.log') or '.cfgt.' in f.path):
        delete_file(f)
