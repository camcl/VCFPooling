#!/bin/bash -l
#SBATCH -A snic2019-8-216
#SBATCH -p node
#SBATCH -N 1
#SBATCH -t 3-00:00:00
#SBATCH -J parallel_pooling
module load python/3.6.8
module load bioinfo-tools
module load bcftools/1.9
module load tabix/0.2.6
source ~/1000Genomes/venv3.6/bin/activate
python3 ~/1000Genomes/scripts/VCFPooling/python/parallel_20200222.py ~/1000Genomes/data/gt/ALL.chr20.snps.gt.vcf.gz ~/1000Genomes/data/gl/gl_adaptive/parallel_pooling/ALL.chr20.pooled.snps.gl.vcf.gz 20

# adaptive GL pooling with EM method values
# whole chromosom 20 data set filtered for SNP-only

# 2020-02-26 11:00
# [camcl609@rackham1 slurm]$ sbatch parallel_pooling_20200226.sh
# Submitted batch job 12220802
