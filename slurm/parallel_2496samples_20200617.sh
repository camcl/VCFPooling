#!/bin/bash -l
#SBATCH -A snic2019-8-216
#SBATCH -p node
#SBATCH -N 1
#SBATCH -t 00:15:00
#SBATCH -J chr20_2496samples
module load python/3.6.8
module load bioinfo-tools
module load bcftools/1.9
module load tabix/0.2.6
echo 'modules loaded'
bash bcf2496para.sh /crex/proj/snic2019-8-216/chr20_chunks /crex/proj/snic2019-8-216/chr20_2496_samples/ /crex/proj/snic2019-8-216/chr20_2496_samples/ALL.chr20.snps.allID.txt

# Job for extracting 2496 samples from each chunk of chr20 VCF-files
# Based on a parallelized bash script using bcftools for processing VCF files
