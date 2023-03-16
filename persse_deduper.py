import argparse
import re
import pysam

parser = argparse.ArgumentParser()
def getArgs():
	parser.add_argument(
		"-f", "--alignment_file", help="Absolute file path to sorted input SAM file", type=str, required=True
	)
	parser.add_argument(
		"-p", "--paired", help="If flag is set, allows for paired end sam deduping if enabled. Paired end functionality is not yet included in this script.", action="store_true", required=False
	)
	parser.add_argument(
		"-u", "--umi", help="Absolute path to file containing UMIs, if provided", type=str, required=False, default=None
		)
	return parser.parse_args()

args = getArgs()

### defining functions ###

def convert_phred(letter: str, ch=33):
	"""Converts a single character into a phred score"""
	return ord(letter) - ch

if __name__ =="main":
	assert convert_phred("I") == 40, "wrong phred score for 'I'"
	assert convert_phred("2") == 17, "wrong phred score for '2'"
	assert convert_phred("@") == 31, "wrong phred score for '@'"

def mean_qual_score(phred_score):
	"""Converts a quality score line into a mean quality score"""
	x = 0
	for ch in phred_score:
		x += (convert_phred(ch))
	return (x / len(phred_score))

if __name__ == "main":
	assert mean_qual_score('GGGGTTTTTGGGGGGG') == 42.0625
	assert mean_qual_score('GOOGOOGAGA') == 40.0
	assert mean_qual_score('!!!!!!!!!!') == 0.0

def softclip_fix(cigar, start):
	"""Adjusts the start position in the file to account for softclipping"""
	new_pos = cigar.split('S')
	if new_pos[0].isdigit():
		new = int(start) - int(new_pos[0])
	else:
		new = start
	return(str(new))

if __name__ == "main":
	assert softclip_fix('10S90M', 10) == 0
	assert softclip_fix('3S97M', 10) == 7
	assert softclip_fix('99M1S', 10) == 10

def sam_split(samline):
	"""Parses through a line from a sam file and splits line for later capture"""
	components = samline.split()
	return(components)

if __name__ == "main":
	assert sam_split('x	V	Z') == ['x', 'V', 'Z']

def rev_strand_pos(cigar, start):
	"""Adjusts start position in file to account for reverse strandedness"""
	caught = re.findall(r"\d+\w", cigar)
	start = int(start)
	if "S" in caught[0]:
		caught[0] = "0"
	for element in caught:
		if "I" in element:
			next
		else:
			start += int(re.sub("[A-Z]", "", element))
	return(start)

if __name__ == "main":
	assert rev_strand_pos('10S90M2I900D90S', 0) == 1080
	assert rev_strand_pos('10S90M2I90D90S', 0) == 270
	assert rev_strand_pos('90M2I90D90S', 0) == 270

## establishing args as variables to be opened
randomer = False
f = args.alignment_file # variable to store input file name
p = args.paired # T/F variable to tell if the data is paired end.
if args.umi: # checks to see if UMIs file is included
	u = args.umi
else:
	randomer=True # if randomer true, we know the umis will only be present in headers. 
out = f.split('.bam')
o = out[-2] + '_deduped.bam'

### global variables ###
unique_dict = {}
umis_set = set()
first = True
reverse,forward,uniqe_dict = {},{},{}

### Reading from input files, writing to output files ###

# create AlignmentFile objects for input and output using Pysam
bam_in = pysam.AlignmentFile(f, 'rb')
bam_out = pysam.AlignmentFile(o, 'wb', template=bam_in)

# if we are not starting with an UMIs file, we have to create the UMI_set on the fly... or so I think.
if randomer:
	umi_set = set()
else:
	with open(u, 'r') as umis:
		# open umis file, read through and add each umi to a set of all umis
		for line in umis:
			line=line.strip()
			umis_set.add(line)

		# add all umis to unique dict as key, whose value is a dictionary
	
		for umi in umis_set:
			unique_dict[umi] = {'reverse' : {}, 'forward' : {}}
	
	# open input sam file, write all lines starting with @ to the outpu file
for bam_line in bam_in:
	line=sam_split(bam_line) # split the line into the components desired, see function description above
	header = line[0]
	umi = header.split(':')[-1]
	if randomer:
		umi_set.add(umi)
	flag = int(line[1])
	start = line[3]
	cigar = line[5]
	qual = mean_qual_score(line[10]) # do I use this value?
	chromo = line[2]
	# need to assign current chromo value before entering loop, if this is the first line
	if first == True:
		current_chromo = chromo
		first = False
		# forgot to add the below section in, adds first read to counter
		if ((flag & 16) == 16):
			start = rev_strand_pos(cigar, start)
			if start not in unique_dict[umi]['reverse'].keys():
						unique_dict[umi]['reverse'][start] = bam_line
		else: 
			if 'S' in cigar:
				start = softclip_fix(cigar, start)
			if start not in unique_dict[umi]['forward'].keys():
				unique_dict[umi]['forward'][start] = bam_line

	else:
		if chromo == current_chromo: # if still in the same chromosome, we continue reading in data
			
			if umi in unique_dict.keys():
				if ((flag & 16) == 16):
					start = rev_strand_pos(cigar, start)
					if start not in unique_dict[umi]['reverse'].keys():
						unique_dict[umi]['reverse'][start] = bam_line
				else: 
					if 'S' in cigar:
						start = softclip_fix(cigar, start)
					if start not in unique_dict[umi]['forward'].keys():
						unique_dict[umi]['forward'][start] = bam_line
					
					# optional section to keep highest quality read encountered per umi.
					# if start in unique_dict.keys():
					# 		# optional section to keep highest quality read encountered per umi.
					# 		if qual > mean_qual_score(unique_dict[start][10]):
					# 			unique_dict[start] = line

		else: # if not on the same chromosome, we start writing
			# the writing part
			for entry in unique_dict.values():
				for direction in entry.values():
					for part in direction.values():
						bam_out.write(part + '\n')
			# after writing contents of the dicitonary, reset dicitonary for next chromosome
			for umi in umis_set:
				unique_dict[umi] = {'reverse' : {}, 'forward' : {}}
			print(current_chromo)
			current_chromo = chromo # setting the current chromosome to the new chromosome encountered on this line
			if umi in unique_dict.keys():
				if ((flag & 16) == 16):
					start = rev_strand_pos(cigar, start)
					if start not in unique_dict[umi]['reverse'].keys():
						unique_dict[umi]['reverse'][start] = bam_line
				else: 
					if 'S' in cigar: # checks for forward read softclipping
						start = softclip_fix(cigar, start)
					if start not in unique_dict[umi]['forward'].keys():
						unique_dict[umi]['forward'][start] = bam_line
# for the last portion which is excluded from our for loop above, since the final chromosome in the sorted file will be excluded from the writing portion.
for entry in unique_dict.values():
	for direction in entry.values():
		for part in direction.values():
			bam_out.write(part + '\n')
