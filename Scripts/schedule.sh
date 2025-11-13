#!/bin/bash
##### Amount of cores per task
#SBATCH --cpus-per-task=32
##### Partition name
#SBATCH -p cpu
##### Name of job in queuing system
#SBATCH --job-name=KTO-DOS
##### name of output file
#SBATCH --output="output.out"
##### name of error file
#SBATCH --error="error.err"

#bash /home/czarnecki/LAO-STO/Scripts/runDoses.sh
srun /home/czarnecki/LAO-STO/bin/LAO_STO.x