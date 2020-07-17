#!/bin/bash -l
#SBATCH -A snic2019-8-216
#SBATCH -p devel
#SBATCH -N 1
#SBATCH --qos=short
#SBATCH -t 00:15:00
#SBATCH -J parallel_pooling
module load python/3.6.8
module load bioinfo-tools
module load bcftools/1.9
module load tabix/0.2.6
echo 'modules loaded'
source ~/1000Genomes/venv3.6/bin/activate
echo 'venv activated'
python3 ~/1000Genomes/scripts/VCFPooling/python/parallel_20200222.py ~/1000Genomes/data/gt/ALL.chr20.snps.gt.chunk10000.vcf.gz ~/1000Genomes/data/gl/gl_adaptive/parallel_pooling/ALL.chr20.pooled.snps.gl.chunk10000.vcf.gz 20

# adaptive GL pooling with EM method values
# 10000 markers stratified per AF-bin
