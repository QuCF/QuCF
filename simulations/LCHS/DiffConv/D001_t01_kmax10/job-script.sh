#!/bin/bash
#SBATCH -N 1
#SBATCH -G 1
#SBATCH -C gpu&hbm80g
#SBATCH -n 1
#SBATCH -c 32
##SBATCH --gpus-per-node=1
#SBATCH --gpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=novikau1@llnl.gov
#SBATCH --account=mp2_g
#SBATCH --qos=shared
#SBATCH -t 00:05:00
#SBATCH -J CDE_t01_kmax10
#
## To run
/global/homes/n/novikau1/QuCF/QuCF/build_qucf/QuCF dc ./
