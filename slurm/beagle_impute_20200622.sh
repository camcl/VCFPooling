#!/bin/bash -l
#SBATCH -A snic2019-8-216
#SBATCH -p devel
#SBATCH --qos=short
#SBATCH -n 20
#SBATCH -t 00:15:00
#SBATCH -J beagle_impute_chr20
module load python/3.6.8
module load bioinfo-tools
module load bcftools/1.9
module load tabix/0.2.6
source ~/1000Genomes/venv3.6/bin/activate
python3 -u /crex/proj/snic2019-8-216/20200622/beagle_impute_20200622.py /crex/proj/snic2019-8-216/20200622 20