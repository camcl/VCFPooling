#!/bin/bash -l

#SBATCH -A snic2019-8-216
#SBATCH -p node
#SBATCH -N 1
#SBATCH -t 04:30:00
#SBATCH -J chr20_pooling
module load python/3.6.8
module load bioinfo-tools
module load bcftools/1.9
module load tabix/0.2.6
echo 'modules loaded'

# pack and write chunks with bcftools and bash
mkdir /crex/proj/snic2019-8-216/chr20_parallel_pooling
tabix /crex/proj/snic2019-8-216/chr20_2496_samples/ALL.chr20.snps.gt.vcf.gz
touch -c /crex/proj/snic2019-8-216/chr20_2496_samples/ALL.chr20.snps.gt.vcf.gz.csi

# pool packed files with python
source ~/1000Genomes/venv3.6/bin/activate
python3 -c 'print("\npython3 ready\n")'
echo ''
python3 -u parallel_pooling_20200617.py /crex/proj/snic2019-8-216/chr20_2496_samples/ALL.chr20.snps.gt.vcf.gz /crex/proj/snic2019-8-216/chr20_parallel_pooling/ALL.chr20.pooled.snps.gl.vcf.gz 20

# adaptive GL pooling with EM method values
# 2496 = 156*16 samples
# whole chromosom 20 data set filtered for SNP-only
