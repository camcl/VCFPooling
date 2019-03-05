import os
import argparse
import subprocess

from scripts.poolSNPs import parameters as prm
from persotools.debugging import *
from persotools.files import *

#TODO: remove miss.pool

'''
Run bcftools manipulations for preprocessing the files which are used for Beagle imputation.
Run Beagle.
'''
chk_sz = prm.CHK_SZ
cd = prm.WD
os.chdir(cd)

raw = prm.RAW

pooled = prm.POOLED

missing = prm.MISSING

miss_pool = prm.MISS_POOL

subset = True #if susbset get the first 1000 lines to recreate ech method file: polled, missing etc
trc = prm.SUBCHUNK

### BGZIP
print('BGZIP in {}'.format(os.getcwd()).ljust(80, '.'))
for dic in [pooled, missing]:
    delete_file(dic['gz'])
    if subset:
        delete_file(dic['gz'].replace('chunk' + str(chk_sz), 'chunk' + str(trc)) + '.vcf.gz')

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

    if subset:
        subprocess.run(bgzip.replace('chunk' + str(chk_sz), 'chunk' + str(trc)), shell=True, cwd=cd)
        subprocess.run(sort.replace('chunk' + str(chk_sz), 'chunk' + str(trc)), shell=True, cwd=cd)
        subprocess.run(idxgz.replace('chunk' + str(chk_sz), 'chunk' + str(trc)), shell=True, cwd=cd)

samples_files = ['cat ALL.chr20.snps.gt.allID.txt | head -{} > ALL.chr20.snps.gt.impID.txt'.format(prm.NB_IMP),
                 'cat ALL.chr20.snps.gt.allID.txt | tail -{} > ALL.chr20.snps.gt.refID.txt'.format(prm.NB_REF),
                 'dos2unix ALL.chr20.snps.gt.refID.txt',
                 'dos2unix ALL.chr20.snps.gt.impID.txt']
for f in samples_files:
    subprocess.run(f, shell=True, cwd=cd)

### REF/IMP SAMPLING
print('REF/IMP SAMPLING'.ljust(80, '.'))
for dic in [raw, pooled, missing]:
    delete_file(dic['imp'])
    if subset:
        delete_file(dic['imp'].replace('chunk' + str(chk_sz), 'chunk' + str(trc)) + '.vcf.gz')

    samp = ' '.join(['bcftools',
                     'view',
                     '-Oz -o',
                     dic['imp'],
                     '-S ALL.chr20.snps.gt.impID.txt',
                     dic['gz']
                     ])

    idxsp = ' '.join(['bcftools',
                      'index -f',
                      dic['imp']
                      ])

    subprocess.run(samp, shell=True, cwd=cd)
    subprocess.run(idxsp, shell=True, cwd=cd)
    if subset:
        subprocess.run(samp.replace('chunk' + str(chk_sz), 'chunk' + str(trc)), shell=True, cwd=cd)
        subprocess.run(idxsp.replace('chunk' + str(chk_sz), 'chunk' + str(trc)), shell=True, cwd=cd)

samp = ' '.join(['bcftools',
                 'view',
                 '-Oz -o',
                 raw['ref'],
                 '-S ALL.chr20.snps.gt.refID.txt',
                 raw['gz']
                 ])

idxsp = ' '.join(['bcftools',
                  '-f index',
                  raw['ref']
                  ])

subprocess.run(samp, shell=True, cwd=cd)
subprocess.run(idxsp, shell=True, cwd=cd)
if subset:
    subprocess.run(samp.replace('chunk' + str(chk_sz), 'chunk' + str(trc)), shell=True, cwd=cd)
    subprocess.run(idxsp.replace('chunk' + str(chk_sz), 'chunk' + str(trc)), shell=True, cwd=cd)

for f in [raw['imp'], raw['ref'], pooled['imp'], missing['imp']]:
    subprocess.run('bcftools sort {} {}'.format(f, f), shell=True, cwd=cd)
    if subset:
        subprocess.run('bcftools sort {} {}'.format(f, f).replace('chunk' + str(chk_sz), 'chunk' + str(trc)),
                       shell=True,
                       cwd=cd)
### BEAGLE ROUND#1
print('BEAGLE ROUND#1'.ljust(80, '.'))
delete_file(raw['b1r'] + '.vcf.gz')
delete_file(raw['b1i'] + '.vcf.gz')
if subset:
    delete_file(raw['b1r'].replace('chunk' + str(chk_sz), 'chunk' + str(trc)) + '.vcf.gz')
    delete_file(raw['b1i'].replace('chunk' + str(chk_sz), 'chunk' + str(trc)) + '.vcf.gz')


bgl1 = ' '.join(['java -Xmx4000m -jar beagle.27Jan18.7e1.jar',
                 'gt=' + raw['imp'],
                 'impute=false',
                 'gprobs=true',
                 'out=' + raw['b1i'],
                 '&',
                 'java -Xmx4000m -jar beagle.27Jan18.7e1.jar',
                 'gt=' + raw['ref'],
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

subprocess.run(bgl1, shell=True, cwd=cd)
subprocess.run(idxb1, shell=True, cwd=cd)
if subset:
    subprocess.run(bgl1.replace('chunk' + str(chk_sz), 'chunk' + str(trc)), shell=True, cwd=cd)
    subprocess.run(idxb1.replace('chunk' + str(chk_sz), 'chunk' + str(trc)), shell=True, cwd=cd)

### CONFORM-GT
#TODO: read pos min et max du fichier chunk
print('CONFORM-GT'.ljust(80, '.'))
for dic in [pooled, missing]:
    delete_file(dic['cfgt'] + '.vcf.gz')
    if subset:
        delete_file(dic['cfgt'].replace('chunk' + str(chk_sz), 'chunk' + str(trc)) + '.vcf.gz')

    cfgt = ' '.join(['java -jar conform-gt.jar',
                     'gt=' + dic['imp'],
                     'chrom=20:60343-62965354',
                     'ref=REF.chr20.snps.gt.chunk{}.vcf.gz'.format(str(chk_sz)),
                     'out=' + dic['cfgt']
                     ])

    idxcf = ' '.join(['bcftools',
                      'index -f',
                      dic['cfgt'] + '.vcf.gz'
                      ])

    subprocess.run(cfgt, shell=True, cwd=cd)
    subprocess.run(idxcf, shell=True, cwd=cd)
    if subset:
        subprocess.run(cfgt.replace('chunk' + str(chk_sz), 'chunk' + str(trc)), shell=True, cwd=cd)
        subprocess.run(idxcf.replace('chunk' + str(chk_sz), 'chunk' + str(trc)), shell=True, cwd=cd)

### BEAGLE (ROUND#2)
print('BEAGLE (ROUND#2)'.ljust(80, '.'))
for dic in [pooled, missing]:
    delete_file(dic['b2'] + '.vcf.gz')
    if subset:
        delete_file(dic['b2'].replace('chunk' + str(chk_sz), 'chunk' + str(trc)) + '.vcf.gz')

    bgl2 = ' '.join(['java -Xmx4000m -jar beagle.27Jan18.7e1.jar',
                     'gt=' + dic['cfgt'] + '.vcf.gz',
                     'ref=REF.chr20.beagle1.chunk{}.vcf.gz'.format(str(chk_sz)),
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
    if subset:
        subprocess.run(bgl1.replace('chunk' + str(chk_sz), 'chunk' + str(trc)), shell=True, cwd=cd)
        subprocess.run(idxb2.replace('chunk' + str(chk_sz), 'chunk' + str(trc)), shell=True, cwd=cd)

### FIX DS AND GP FORMAT FIELDS
print('REFORMATTING GP AND DS FIELDS'.ljust(80, '.'))
for dic in [pooled, missing]:
    delete_file(dic['corr'] + '.vcf.gz')
    if subset:
        delete_file(dic['corr'].replace('chunk' + str(chk_sz), 'chunk' + str(trc)) + '.vcf.gz')

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
    if subset:
        subprocess.run(refmt.replace('chunk' + str(chk_sz), 'chunk' + str(trc)), shell=True, cwd=cd)
        subprocess.run(idxfm.replace('chunk' + str(chk_sz), 'chunk' + str(trc)), shell=True, cwd=cd)


### RESULTS COMPARISON
print('RESULTS plot-vcfstats'.ljust(80, '.'))
for dic in [pooled, missing]:
    straw = ' '.join(['bcftools stats --verbose',
                      'IMP.chr20.beagle1.vcf.gz -t "Raw"',
                      dic['b2'] + '.vcf.gz -t "Processed"',
                      '> ' + dic['b2'] + '.check'
                      ])
    plot = ' '.join(['plot-vcfstats -s ',
                     '-t "Raw vs. {}"'.format(dic['imp']),
                     '-p ./' + dic['b2'],
                     dic['b2'] + '.check'])

    # subprocess.run(straw, shell=True, cwd=cd)
    # subprocess.run(plot, shell=True, cwd=cd)

missdiff = ' '.join(['bcftools stats --verbose',
                     'IMP.chr20.pooled.beagle2.vcf.gz',
                     'IMP.chr20.missing.pooled.beagle2.vcf.gz',
                     '> IMP.chr20.missdiff.check'
                     ])
pooldiff = ' '.join(['bcftools stats --verbose',
                     'IMP.chr20.missing.beagle2.vcf.gz',
                     'IMP.chr20.missing.pooled.beagle2.vcf.gz',
                     '> IMP.chr20.pooldiff.check'
                     ])
pltdiff = ' '.join(['plot-vcfstats -s',
                    '-t "Pooldiff"',
                    '-p ./pooldiff',
                    'IMP.chr20.pooldiff.check',
                    '&',
                    'plot-vcfstats -s',
                    '-t "Missdiff"',
                    '-p ./missdiff',
                    'IMP.chr20.missdiff.check'
                    ])

# subprocess.run(missdiff, shell=True, cwd=cd)
# subprocess.run(pooldiff, shell=True, cwd=cd)
# subprocess.run(pltdiff, shell=True, cwd=cd)