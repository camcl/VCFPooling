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
python3 ~/1000Genomes/scripts/VCFPooling/python/parallel_pooling_20200615.py /crex/proj/snic2019-8-216/20200615/ALL.chr20.snps.gt.chunk10000.vcf.gz /crex/proj/snic2019-8-216/20200615/ALL.chr20.pooled.snps.gl.chunk10000.vcf.gz 20

# adaptive GL pooling with EM method values
# 10000 markers stratified per AF-bin
# new pooling simulation implementation for comparison purposes with 20200222 experiments
