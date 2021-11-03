#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --time=05:00:00
#SBATCH --output=conversion%j.out
#SBATCH --error=conversion%j.err
#SBATCH --mail-user='tpersse@uoregon.edu'
#SBATCH --mail-type=END,FAIL

/usr/bin/time -v samtools sort C1_SE_uniqAlign.sam -o sorted_C1_SE.sam