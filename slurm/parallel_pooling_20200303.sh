#!/bin/bash -l

#SBATCH -A snic2019-8-216
#SBATCH -p node
#SBATCH -N 1
#SBATCH -t 3-00:00:00
#SBATCH -J chr20_pooling
module load python/3.6.8
module load bioinfo-tools
module load bcftools/1.9
module load tabix/0.2.6
echo 'modules loaded'

# pack and write chunks with bcftools and bash
mkdir /proj/snic2019-8-216/parallel_vcfsplit
cd /proj/snic2019-8-216/parallel_vcfsplit
cp ~/1000Genomes/data/gt/ALL.chr20.snps.gt.vcf.gz /proj/snic2019-8-216/parallel_vcfsplit
tabix ALL.chr20.snps.gt.vcf.gz
touch -c ALL.chr20.snps.gt.vcf.gz.csi
bash ~/1000Genomes/scripts/bcfchunkpara.sh ALL.chr20.snps.gt.vcf.gz

# pool packed files with python
source ~/1000Genomes/venv3.6/bin/activate
python3 -c 'print("\npython3 ready\n")'
echo ''
python3 -u ~/1000Genomes/scripts/VCFPooling/python/parallel_20200303.py /crex/proj/snic2019-8-216/parallel_vcfsplit/ALL.chr20.snps.gt.vcf.gz ~/1000Genomes/data/gl/gl_adaptive/parallel_pooling/ALL.chr20.pooled.snps.gl.vcf.gz

# adaptive GL pooling with EM method values
# whole chromosom 20 data set filtered for SNP-only
