#!/bin/bash -l
#SBATCH -A snic2019-8-216
#SBATCH -p node
#SBATCH -N 1
#SBATCH -t 03:00:00
#SBATCH -J chr20_splitting
module load python/3.6.8
module load bioinfo-tools
module load bcftools/1.9
module load tabix/0.2.6
echo 'modules loaded'
mkdir /proj/snic2019-8-216/parallel_vcfsplit
cd /proj/snic2019-8-216/parallel_vcfsplit
cp ~/1000Genomes/data/gt/ALL.chr20.snps.gt.vcf.gz /proj/snic2019-8-216/parallel_vcfsplit
tabix ALL.chr20.snps.gt.vcf.gz
touch -c ALL.chr20.snps.gt.vcf.gz.csi
bash ~/1000Genomes/src/bcfchunkpara.sh ALL.chr20.snps.gt.vcf.gz

# Job for splitting a chr20 VCF-file into chunks of 1000 markers
# Based on a parallelized bash script using bcftools for processing VCF files
