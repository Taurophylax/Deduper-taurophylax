#!/bin/bash
#SBATCH --partition=bgmp       ### Partition (like a queue in PBS)
#SBATCH --time=0-12:00:00       ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1               ### Number of nodes needed for the job
#SBATCH --ntasks-per-node=1     ### Number of tasks to be launched per Node
#SBATCH --cpus-per-task=8
#SBATCH --account=bgmp         ### Account used for job submission
#SBATCH --mem=28G

/usr/bin/time -v python taurophylax_deduper.py -u STL96.txt -f input.sam -o output.sam
