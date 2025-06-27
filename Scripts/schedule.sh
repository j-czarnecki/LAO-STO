#!/bin/bash
##### Amount of cores per task
#SBATCH --cpus-per-task=64
##### Partition name
#SBATCH -p cpu
##### Name of job in queuing system
#SBATCH --job-name=LAO-STO

#bash /home/czarnecki/LAO-STO/Scripts/runDoses.sh
srun -c 48 /home/czarnecki/LAO-STO/bin/POST_LAO_STO.x
make analyze