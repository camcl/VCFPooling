import subprocess
import os
import pathlib
from cyvcf2 import VCF

from scripts.VCFPooling.poolSNPs import parameters as prm
from scripts.VCFPooling.poolSNPs.myssh import sshtools
from scripts.VCFPooling.poolSNPs import pybcf
from persotools.files import delete_file

'''
Phaser C++ code run on pooled samples for 10000 stratified markers.
Imputation on GL filled with adaptative values.
1 file per sample as output, with 1 relevant column
Files location:
/proj/snic2019-8-164/private/Phaser/Camille_out6

1) Select the relevant column for each sample: bcftools view -Oz -o <file_out> -s<ID_SAMPLE> <file_in> 
2) Merge the columns in 1 vcf.gz file named: prm.POOLED['gt-only']
3) Copy the remote directory to local directory gl/gl_adaptative/phaser
4) Plot results
'''

remote = False
local = False
reorder = True

# Paths/Filenames settings
localpath = pathlib.PurePath('/home/camille/1000Genomes')
littlelocal = pathlib.PurePath('/home/camille/1000G/')
remotepath = pathlib.PurePath('/proj/snic2019-8-164/private/Phaser/Camille_out_190822')
remote_fname = 'IMP.chr20.pooled.snps.gl.chunk10000.full_{}.genos.vcf.gz'
# Read and index samples IDs from the file send
inputpath = pathlib.PurePath('/home/camille/1000Genomes/data/chunk10000_stratifiedAAF_20190822')
inputfile = 'IMP.chr20.pooled.snps.gl.chunk10000.vcf.gz'
outputpath = pathlib.PurePath(localpath, 'data', 'gl', 'gl_adaptative', 'phaser')
sampleslist = VCF(str(pathlib.PurePath(inputpath, inputfile))).samples
sampleslist = list([s.strip('_IMP') for s in sampleslist])
# e.g. HG01125_IMP to HG01125, no suffix in remote folder
samplestab = set(zip(sampleslist, range(len(sampleslist))))


### Connect to the remote repository
if remote:
    with sshtools.ssh_connect(host=sshtools.host,
                              port=22,
                              user=sshtools.user,
                              connect_kwargs={'password': sshtools.pwd}) as client:
        print('\nRemote processing on {}:{}'.format(sshtools.host, remotepath))
        ls1 = client.run('ls', pty=True)
        cdls1 = client.run('cd {} && ls'.format(remotepath), pty=True)
        load_bcftls = 'module load bioinfo-tools && module load bcftools/1.8'
        files = list()
        for spl, idx in samplestab:
            cmd = ' '.join(['bcftools',
                            'view -Oz -o',
                            'tmp.{}'.format(remote_fname.format(idx)),
                            '-s {}'.format(spl),
                            remote_fname.format(idx),
                            '&&',
                            'bcftools index -f tmp.{}'.format(remote_fname.format(idx))
                            ])
            extract = client.run('cd {} && {} && {}'.format(remotepath, load_bcftls, cmd),
                                 pty=True)
            files.append('tmp.{}'.format(remote_fname.format(idx)))
            print(cmd)
        files_str = ' '.join(files)
        print('Merge files: {}'.format(files_str))
        bcfmerge = ' '.join(['bcftools merge',
                             '-Oz -o',
                             prm.POOLED['gtonly'] + '.vcf.gz',
                             files_str
                             ])
        merge = client.run('cd {} && {} && {}'.format(remotepath, load_bcftls, bcfmerge),
                           pty=True)
        clean = client.run('cd {} && rm tmp.*'.format(remotepath), pty=True)
        cdls2 = client.run('cd {} && ls'.format(remotepath), pty=True)

# Close connection! Run LOCALLY:
if local:
    print('\nLocal processing in {}'.format(localpath))
    mkphaser = ' '.join(['cd',
                         os.path.join(str(localpath), 'data', 'gl', 'gl_adaptative'),
                         '&&'])
                         # 'mkdir phaser'])
    scp = ' '.join(['scp',
                    '{}:{}'.format(sshtools.uppmax, str(pathlib.PurePath(remotepath, prm.POOLED['gtonly'])) + '.vcf.gz'),
                    os.path.join(str(outputpath), prm.POOLED['gtonly'] + '.vcf.gz')])
    subprocess.run(mkphaser, shell=True)
    subprocess.run(scp, shell=True)
    # subprocess.run(' '.join(['cp -r',
    #                          str(inputpath) + '/*.vcf.gz*',
    #                          str(outputpath) + '/*.vcf.gz*']))

    # Index the file
    pybcf.index(prm.POOLED['gtonly'] + '.vcf.gz',
                str(outputpath))

# # Recreate GT pooled from GL pooled: 2 times
# bgl1gt = ' '.join(['java -Xss5m -jar {}'.format(prm.BEAGLE_JAR),
#                    '{}='.format('gtgl') + prm.POOLED['imp'],
#                    'impute=false',
#                    'gprobs=true',
#                    'out=' + prm.POOLED['imp'].strip('.vcf.gz').replace('.gl', '.gt_unphased')])
# bgl2gt = ' '.join(['java -Xss5m -jar {}'.format(prm.BEAGLE_JAR),
#                    '{}='.format('gt') + prm.POOLED['imp'].replace('.gl', '.gt_unphased'),
#                    'impute=false',
#                    'gprobs=true',
#                    'out=' + prm.POOLED['imp'].replace('.gl', '.gt').strip('.vcf.gz')
#                    ])
# print(bgl1gt)
# print(bgl2gt)
# subprocess.run(bgl1gt, shell=True, cwd=str(outputpath))
# subprocess.run(bgl2gt, shell=True, cwd=str(outputpath))
# pybcf.index(prm.POOLED['imp'].replace('.gl', '.gt'), str(outputpath))
#

# Reorder samples according to pools
if reorder:
    print('\nReorder imputed samples locally in {}'.format(outputpath))
    os.chdir(outputpath)
    subprocess.run(' '.join(['bcftools view',
                             '-Oz -o IMP.chr20.pooled.imputed.chunk{}.unordered.samples.vcf.gz'.format(str(prm.CHK_SZ)),
                             'IMP.chr20.pooled.imputed.gt.chunk{}.vcf.gz'.format(str(prm.CHK_SZ))]),
                   shell=True,
                   cwd=outputpath)

    subprocess.run('bcftools index -f IMP.chr20.pooled.imputed.chunk{}.unordered.samples.vcf.gz'.format(str(prm.CHK_SZ)),
                   shell=True,
                   cwd=outputpath)

    subprocess.run(' '.join(['bcftools view -S {}/ALL.chr20.snps.impID.txt'.format(prm.WD + '/gt'),
                             '-Oz -o IMP.chr20.pooled.imputed.gt.chunk{}.vcf.gz'.format(str(prm.CHK_SZ)),
                             'IMP.chr20.pooled.imputed.chunk{}.unordered.samples.vcf.gz'.format(str(prm.CHK_SZ))]),
                   shell=True,
                   cwd=outputpath)

    subprocess.run('bcftools index -f IMP.chr20.pooled.imputed.gt.chunk{}.vcf.gz'.format(str(prm.CHK_SZ)),
                   shell=True,
                   cwd=outputpath)

    delete_file('IMP.chr20.pooled.imputed.chunk{}.unordered.samples.vcf.gz'.format(str(prm.CHK_SZ)))
    delete_file('IMP.chr20.pooled.imputed.chunk{}.unordered.samples.vcf.gz.csi'.format(str(prm.CHK_SZ)))
