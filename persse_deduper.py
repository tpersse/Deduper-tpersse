import argparse
import re


parser = argparse.ArgumentParser()
def getArgs():
	parser.add_argument(
		"-f", "--file", help="Absolute file path to sorted input SAM file", type=str, required=True
	)
	parser.add_argument(
		"-p", "--paired", help="If flag is set, allows for paired end sam deduping if enabled. Paired end functionality is not yet included in this script.", action="store_true", required=False
	)
	parser.add_argument(
		"-u", "--umi", help="Absolute path to file containing UMIs, if provided", action="store_true", type=str, required=True
	)

	return parser.parse_args()

args = getArgs()
f = args.file
p = args.paired
u = args.umi


if p:
	parser.error("paired end functionality has not yet been built into this script")


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
	assert qual_score('GGGGTTTTTGGGGGGG') == 42.0625
	assert qual_score('GOOGOOGAGA') == 40.0
	assert qual_score('!!!!!!!!!!') == 0.0

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
	components = samline.strip().split()
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




### creating output file name ###

out = f.split('.sam')
o = out[-2] + '_deduped.sam'

### global variables ###

unique_dict = {}

umis_set = set()

first = True

reverse = {}

forward = {}


### Reading from input files, writing to output files ###

with open(u, 'r') as umis, open(f, 'r') as sam_in, open(o, 'w') as sam_out:
	# open umis file, read through and add each umi to a set of all umis
	for line in umis:
		line=line.strip()
		umis_set.add(line)

	# add all umis to unique dict as key, whose value is a dictionary
	for umi in umis_set:
		unique_dict[umi] = {'reverse' : {}, 'forward' : {}}
	
	# open input sam file, write all lines starting with @ to the outpu file
	for sam_line in sam_in:
		sam_line = sam_line.strip()
		if sam_line.startswith('@'):
			sam_out.write(sam_line + '\n')
		else: # any line in sam file that contains sequence
			line=sam_split(sam_line) # split the line into the components desired, see function description above
			header = line[0]
			umi = header.split(':')[-1]
			flag = int(line[1])
			start = line[3]
			cigar = line[5]
			qual = mean_qual_score(line[10])
			chromo = line[2]
			# need to assign current chromo value before entering loop, if this is the first line
			if first == True:
				current_chromo = chromo
				first = False
				# forgot to add the below section in, adds first read to counter
				if ((flag & 16) == 16):
					start = rev_strand_pos(cigar, start)
					if start not in unique_dict[umi]['reverse'].keys():
								unique_dict[umi]['reverse'][start] = sam_line
				else: 
					if 'S' in cigar:
						start = softclip_fix(cigar, start)
					if start not in unique_dict[umi]['forward'].keys():
						unique_dict[umi]['forward'][start] = sam_line

			else:
				if chromo == current_chromo: # if still in the same chromosome
					
					if umi in unique_dict.keys():
						if ((flag & 16) == 16):
							start = rev_strand_pos(cigar, start)
							if start not in unique_dict[umi]['reverse'].keys():
								unique_dict[umi]['reverse'][start] = sam_line
						else: 
							if 'S' in cigar:
								start = softclip_fix(cigar, start)
							if start not in unique_dict[umi]['forward'].keys():
								unique_dict[umi]['forward'][start] = sam_line
							# optional section to keep highest quality read encountered per umi.
							# if start in unique_dict.keys():
							# 		# optional section to keep highest quality read encountered per umi.
							# 		if qual > mean_qual_score(unique_dict[start][10]):
							# 			unique_dict[start] = line
				
				else: # if not on the same chromosome
					# the writing part
					for entry in unique_dict.values():
						for direction in entry.values():
							for part in direction.values():
								sam_out.write(part + '\n')
					# after writing contents of the dicitonary, reset dicitonary for next chromosome
					for umi in umis_set:
						unique_dict[umi] = {'reverse' : {}, 'forward' : {}}
					print(current_chromo)
					current_chromo = chromo # setting the current chromosome to the new chromosome encountered on this line
					if umi in unique_dict.keys():
						if ((flag & 16) == 16):
							start = rev_strand_pos(cigar, start)
							if start not in unique_dict[umi]['reverse'].keys():
								unique_dict[umi]['reverse'][start] = sam_line
						else: 
							if 'S' in cigar: # checks for forward read softclipping
								start = softclip_fix(cigar, start)
							if start not in unique_dict[umi]['forward'].keys():
								unique_dict[umi]['forward'][start] = sam_line
	for entry in unique_dict.values():
		for direction in entry.values():
			for part in direction.values():
				sam_out.write(part + '\n')