import os
import argparse
import subprocess

'''
Run bcftools manipulations for preprocessing the files which are used for Beagle imputation.
Run Beagle.
'''

os.chdir('/home/camilleclouard/PycharmProjects/1000Genomes/data/tests-beagle')
cd = '/home/camilleclouard/PycharmProjects/1000Genomes/data/tests-beagle'

raw = {'vcf':'ALL.chr20.snps.gt.chunk.vcf',
       'gz':'ALL.chr20.snps.gt.chunk.vcf.gz',
       'ref': 'REF.chr20.snps.gt.chunk.vcf.gz',
        'imp': 'IMP.chr20.snps.gt.chunk.vcf.gz',
       'b1r':'REF.chr20.beagle1',
       'b1i':'IMP.chr20.beagle1'}

pooled = {'vcf':'ALL.chr20.pooled.snps.gt.chunk.vcf',
       'gz':'ALL.chr20.pooled.snps.gt.chunk.vcf.gz',
       'imp': 'IMP.chr20.pooled.snps.gt.chunk.vcf.gz',
       'b1':'IMP.chr20.pooled.beagle1',
       'b2':'IMP.chr20.pooled.beagle2',
          'corr':'IMP.chr20.pooled.beagle2.corr',
        'cfgt': 'IMP.chr20.pooled.cfgt'}

missing = {'vcf':'ALL.chr20.missing.snps.gt.chunk.vcf',
       'gz':'ALL.chr20.missing.snps.gt.chunk.vcf.gz',
       'imp': 'IMP.chr20.missing.snps.gt.chunk.vcf.gz',
       'b1':'IMP.chr20.missing.beagle1',
       'b2':'IMP.chr20.missing.beagle2',
           'corr': 'IMP.chr20.missing.beagle2.corr',
           'cfgt': 'IMP.chr20.missing.cfgt'}

miss_pool = {'vcf':'ALL.chr20.missing.snps.gt.chunk.vcf',
       'gz':'ALL.chr20.missing.pooled.snps.gt.chunk.vcf.gz',
       'imp': 'IMP.chr20.missing.pooled.snps.gt.chunk.vcf.gz',
       'b1':'IMP.chr20.missing.pooled.beagle1',
       'b2':'IMP.chr20.missing.pooled.beagle2',
             'corr': 'IMP.chr20.missing.pooled.beagle2.corr',
             'cfgt': 'IMP.chr20.missing.pooled.cfgt'}


def delete_file(file_path):
    """
    Deletes an existing local file
    :param file_path: string
    :return: -
    """
    if os.path.exists(file_path):
        os.remove(file_path)
    else:
        print("The file does not exists")


### BGZIP
print('BGZIP in {}'.format(os.getcwd()).ljust(80, '.'))
for dic in [pooled, missing, miss_pool]:
    delete_file(dic['gz'])

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

samples_files = ['bcftools query -l ALL.chr20.snps.gt.chunk.vcf.gz | shuf > ALL.chr20.snps.gt.allID.txt',
                'cat ALL.chr20.snps.gt.allID.txt | head -250 > ALL.chr20.snps.gt.impID.txt',
                'cat ALL.chr20.snps.gt.allID.txt | tail -2246 > ALL.chr20.snps.gt.refID.txt',
                'dos2unix ALL.chr20.snps.gt.refID.txt',
                'dos2unix ALL.chr20.snps.gt.impID.txt']
for f in samples_files:
    subprocess.run(f, shell=True, cwd=cd)

### REF/IMP SAMPLING
# TODO: problem when calling sampling from Python, probably because of bash location
print('REF/IMP SAMPLING'.ljust(80, '.'))
for dic in [raw, pooled, missing, miss_pool]:
    delete_file(dic['imp'])

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

for f in [raw['imp'], raw['ref'], pooled['imp'], missing['imp'], miss_pool['imp']]:
    subprocess.run('bcftools sort {} {}'.format(f, f), shell=True, cwd=cd)

### BEAGLE ROUND#1
print('BEAGLE ROUND#1'.ljust(80, '.'))
delete_file(raw['b1r'] + '.vcf.gz')
delete_file(raw['b1i'] + '.vcf.gz')

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

### CONFORM-GT
#TODO: read pos min et max du fichier chunk
print('CONFORM-GT'.ljust(80, '.'))
for dic in [pooled, missing, miss_pool]:
    delete_file(dic['cfgt'] + '.vcf.gz')

    cfgt = ' '.join(['java -jar conform-gt.jar',
                    'gt=' + dic['imp'],
                     'chrom=20:60343-62965354',
                     'ref=REF.chr20.snps.gt.chunk.vcf.gz',
                    'out=' + dic['cfgt']
                    ])

    idxcf = ' '.join(['bcftools',
                      'index -f',
                      dic['cfgt'] + '.vcf.gz'
                      ])

    subprocess.run(cfgt, shell=True, cwd=cd)
    subprocess.run(idxcf, shell=True, cwd=cd)

### BEAGLE (ROUND#2)
print('BEAGLE (ROUND#2)'.ljust(80, '.'))
for dic in [pooled, missing, miss_pool]:
    delete_file(dic['b2'] + '.vcf.gz')

    bgl2 = ' '.join(['java -Xmx4000m -jar beagle.27Jan18.7e1.jar',
                     'gt=' + dic['cfgt'] + '.vcf.gz',
                     'ref=REF.chr20.beagle1.vcf.gz',
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

### FIX DS AND GP FORMAT FIELDS
print('REFORMATTING GP AND DS FIELDS'.ljust(80, '.'))
for dic in [pooled, missing, miss_pool]:
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



### RESULTS COMPARISON
print('RESULTS plot-vcfstats'.ljust(80, '.'))
for dic in [pooled, missing, miss_pool]:
    straw = ' '.join(['bcftools stats --verbose',
                     'IMP.chr20.beagle1.vcf.gz -t "Raw"',
                     dic['b2'] + '.vcf.gz -t "Processed"',
                      '> ' + dic['b2'] + '.check'
                     ])
    plot = ' '.join(['plot-vcfstats -s ',
                     '-t "Raw vs. {}"'.format(dic['imp']),
                     '-p ./' + dic['b2'],
                     dic['b2'] + '.check'])

    subprocess.run([straw, plot], shell=True, cwd=cd)

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

subprocess.run([missdiff, pooldiff, pltdiff], shell=True, cwd=cd)
