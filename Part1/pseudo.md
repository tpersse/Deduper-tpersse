# Assignment the First:

## The Problem:
Within a given fastq/SAM file (any file which contains sequencing reads that have been amplified via PCR), there exist duplicate reads due to library amplification. 
In cases like transcriptomics and genome assembly, this is a vital step, since identical reads could give the impression of increased expression of a gene in the case of the former scenario, or repeat regions in the case of the latter.

Input files: 
 SAM file with uniquely mapped reads, list of UMIs

From slides: 

What does PCR-duplicate look like?

Same alignment position 

•Chromosome 

•Position 

•Strand (strand specific?) 

•Soft Clipping 

•Same Unique Molecular Index (UMI or “randomer”) 

•Single-end vs Paired-end?

general structure of sorting algorithm:

### Global Variable Establishment, List Building, etc.

chromos_dict: a dictionary whose keys will be the chromosome numbers present in the data, and whose values will open files respective to each chromosome [one for writing to, one for reading from] (a with open statement)

### Reading in UMIs
read in UMIs:
create list of UMIs, for indexing later

### SAM file sorting
using for loop, read through sam file line by line:
	split line via split()
	if column 3 (chromosome number) is in chromos_dict.keys():
		write the line to the file
	else:
		create entry in chromos_dict, where key is column 3 for this line, and values are a file to write to (position 0), and a file to read from(position 1)

open an output file for writing to, able to write to it mutliple times (not with 'w' option)

	create a for loop, loop through chromos_dict.values(), position 1 (opens the files in read mode):

		create empty dicitonary, UMIs_dict. Will reset every time a new file is opened, keys are UMIs, values are list of lines with that UMI

		loop through file line by line using a for loop: 

			split line

			if the UMI is in UMIs_dict(position 0 in split line, but the last part of the portion that was split, use regex to assess)
				
				if read on same strand (assessed by looking at bitwise flag):

					if same start AND stop position:
						PCR duplicate, so keep strand with best quality score sequence as dictionary value for UMI (position 10)(acheived via comparison of average phred33 score, function in bioinfo.py)

					else:
						if start position != the same:

							index cigar string and determine where softclipping occurred, and...

								view cigar string of line, determine how many entries are softclipped(x nt clipped)
								
								use string slicing to compare current line to dictionary entry
								
								if line != dictionary value after this comparison:
									
									write the current line to the output sam file (it has passed the sorting algorithm, is softclipped, but is not a duplicate and can be added to the output file)
				
				else:
					
					take reverse complement of sequence, then run it through last part of sorting algorithm:
					
					if same start AND stop position:
						
						PCR duplicate, so keep strand with best quality score sequence as dictionary value for UMI (position 10)(acheived via comparison of average phred33 score, function in bioinfo.py)
					
					else:
						
						if start position != the same:
							
							index cigar string and determine where softclipping occurred, and...
								
								view cigar string of line, determine how many entries are softclipped(x nt clipped)
								
								use string slicing to compare current line to dictionary entry
								
								if line != dictionary value after this comparison:
								
									write the current line to the output sam file (it has passed the sorting algorithm, is softclipped, but is not a duplicate and can be added to the output file)
			
			else:
				
				if UMI in list of UMIs:
					
					generate entry in UMIs_dict, where UMI is key and value is the line from the SAM file
				
				else:
					
					throw it out
		
		at end of file, write UMI_dict values to output SAM file


## High Level Functions:

def reverse_comp():

```This function takes a sequence of DNA bases and converts it to the reverse complement of the sequence```

	return(rev_comp)

ex. 

in: reverse_comp('ATGTC')

out: 'GACAT'


quality score average:

def avg_qual_score():

```This function takes a string of ASCII characters and returns the average of their phred scores```

	return(avg_qual)

ex. 

in: avg_qual_score('AAAIII')

out: 36