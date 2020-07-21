#!/bin/bash -l
#SBATCH -A snic2019-8-216
#SBATCH -p devel
#SBATCH --qos=short
#SBATCH -n 20
#SBATCH -t 00:15:00
#SBATCH -J beagle_impute_v1.3
module load python/3.6.8
module load bioinfo-tools
module load bcftools/1.9
module load tabix/0.2.6
source ~/1000Genomes/venv3.6/bin/activate
python3 -u ~/1000Genomes/scripts/VCFPooling/python/beagle_impute_v1.3.py chr20_20200307