#!/bin/bash
##### Amount of cores per task
#SBATCH --cpus-per-task=32
##### Partition name
#SBATCH -p cpu
##### Name of job in queuing system
#SBATCH --job-name=KTO-DOS

#bash /home/czarnecki/LAO-STO/Scripts/runDoses.sh
srun /home/czarnecki/LAO-STO/bin/POST_LAO_STO.x