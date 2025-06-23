#!/bin/bash
##### Amount of cores per task
#SBATCH --cpus-per-task=16
##### Partition name
#SBATCH -p cpu
##### Name of job in queuing system
#SBATCH --job-name=LAO-STO

srun bin/LAO_STO.x