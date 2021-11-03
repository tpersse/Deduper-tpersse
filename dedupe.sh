#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --time=05:00:00
#SBATCH --output=dedupe%j.out
#SBATCH --error=dedupe%j.err
#SBATCH --mail-user='tpersse@uoregon.edu'
#SBATCH --mail-type=END,FAIL

/usr/bin/time -v python persse_deduper.py  -f ./sorted_C1_SE.sam -u ./STL96.txt