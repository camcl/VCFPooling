import subprocess

from scripts.poolSNPs import parameters as prm
from scripts.poolSNPs.alleles import alleles_tools as alltls
from persotools.debugging import *
from persotools.files import *

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

subset = prm.SUBSET #if susbset get the first 1000 lines to recreate ech method file: pooled, missing etc
trc = prm.SUBCHUNK
GTGL = 'GT'

### GL CONVERSION
if GTGL == 'GL':
    print('GT to GL in {}'.format(os.getcwd()).ljust(80, '.'))
    for dic in [pooled, missing, raw]:
        alltls.file_likelihood_converter(dic['vcf'], dic['vcf'].replace('.gt', '.gl'))
        for k, v in dic.items():
            dic[k] = v.replace('gt', 'gl')

### BGZIP
print('BGZIP in {}'.format(os.getcwd()).ljust(80, '.'))
for dic in [pooled, missing]:
    delete_file(dic['gz'])
    delete_file(dic['gz'] + '.csi')
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

samples_files = ['cat ALL.chr20.snps.allID.txt | head -{} > ALL.chr20.snps.impID.txt'.format(prm.NB_IMP),
                 'cat ALL.chr20.snps.allID.txt | tail -{} > ALL.chr20.snps.refID.txt'.format(prm.NB_REF),
                 'dos2unix ALL.chr20.snps.refID.txt',
                 'dos2unix ALL.chr20.snps.impID.txt']
for f in samples_files:
    subprocess.run(f, shell=True, cwd=cd)

### REF/IMP SAMPLING
print('REF/IMP SAMPLING'.ljust(80, '.'))
for dic in [raw, pooled, missing]:
    delete_file(dic['imp'])
    delete_file(dic['imp'] + '.csi')
    if subset:
        delete_file(dic['imp'].replace('chunk' + str(chk_sz), 'chunk' + str(trc)) + '.vcf.gz')

    samp = ' '.join(['bcftools',
                     'view',
                     '-Oz -o',
                     dic['imp'],
                     '-S ALL.chr20.snps.impID.txt',
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
                 '-S ALL.chr20.snps.refID.txt',
                 raw['gz']
                 ])

idxsp = ' '.join(['bcftools',
                  '-f index',
                  raw['ref']
                  ])

delete_file(raw['ref'] + '.csi')
subprocess.run(samp, shell=True, cwd=cd)
subprocess.run(idxsp, shell=True, cwd=cd)
if subset:
    delete_file(raw['ref'].replace('chunk' + str(chk_sz), 'chunk' + str(trc)) + '.csi')
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


bgl1 = ' '.join(['java -Xmx4000m -jar {}'.format(prm.BEAGLE_JAR),
                 '{}='.format('gt' if GTGL == 'GT' else 'gl') + raw['imp'],
                 'impute=false',
                 'gprobs=true',
                 'out=' + raw['b1i'],
                 '&',
                 'java -Xmx4000m -jar {}'.format(prm.BEAGLE_JAR),
                 '{}='.format('gt' if GTGL == 'GT' else 'gl') + raw['ref'],
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

    cfgt = ' '.join(['java -jar {}'.format(prm.CFGT_JAR),
                     '{}='.format('gt' if GTGL == 'GT' else 'gl') + dic['imp'],
                     'chrom=20:60343-62965354',
                     'ref={}'.format(str(raw['ref'])),
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

    bgl2 = ' '.join(['java -Xmx4000m -jar {}'.format(prm.BEAGLE_JAR),
                     '{}='.format('gt' if GTGL == 'GT' else 'gl') + dic['cfgt'] + '.vcf.gz',
                     'ref={}.vcf.gz'.format(str(raw['b1r'])),
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
        subprocess.run(bgl2.replace('chunk' + str(chk_sz), 'chunk' + str(trc)), shell=True, cwd=cd)
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
