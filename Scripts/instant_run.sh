#!/bin/bash -l
## Job name
#SBATCH -J Fit-J1-30mT-fine-tuning
## Number of allocated nodes
#SBATCH -N 1
## Number of tasks per node (by default this corresponds to the number of cores allocated per node)
#SBATCH --ntasks-per-node=48
## Memory allocated per core (default is 5GB), comment if mem for whole job should be taken
##SBATCH --mem-per-cpu=3800MB
## Memory allocated for whole job, comment if mem-per-cpu should be taken
#SBATCH --mem=20GB
## Max task execution time (format is HH:MM:SS)
#SBATCH --time=72:00:00
## Name of grant to which resource usage will be charged
#SBATCH -A plglaosto111-cpu
## Name of partition
#SBATCH -p plgrid
## Name of file to which standard output will be redirected
#SBATCH --output="output_J1_30mT.out"
## Name of file to which the standard error stream will be redirected
#SBATCH --error="error_J1_30mT.err"

#cd Analyzer
#python3 mainAnalyzer.py
python3 Analyzer/DosFitter.py
