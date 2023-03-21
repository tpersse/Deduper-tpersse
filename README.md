# Deduper
### Author: Thomas Persse

---

Script for the removal of PCR duplicates from an alignment file file. Returns the highest quality read (based on average sequence quality score) for every unique UMI + start position combination. 

Can be given an UMI file for reference if operating with a known pool of UMIs, or can operate UMI-unaware. 

Paired end functionality is still a work in progress. 

## Usage:
` perse_deduper.py [-f --alignment_file] [-u --umi_file] [-p --paired_end]`

## Dependencies:
python=3.9
pysam=0.20.0
All dependencies can be found in the attached .yml file. 
Note: .yml file contains additional package, samtools, for use in the preprocessing of data entering the pipeline.
