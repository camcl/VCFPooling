#!/bin/bash -l
#SBATCH -A snic2019-8-216
#SBATCH -p devel
#SBATCH --qos=short
#SBATCH -n 4
#SBATCH -t 00:02:00
#SBATCH -J hello_world
module load python/3.6.8
module load bioinfo-tools
module load bcftools/1.9
module load tabix/0.2.6
source ~/1000Genomes/venv3.6/bin/activate
module load openmpi/4.0.2
pip install mpi4py
mpirun -n 4 python3 -u ~/PDC-MPI/1_Intro/hello.py