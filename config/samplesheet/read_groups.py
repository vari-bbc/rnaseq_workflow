import os
import gzip
import subprocess

# The goal of this script is to build read groups into units.tsv IF NEEDED 

# Those files are in the raw_read files

# Set up some variables that can be easily changed later
path_to_raw_data = "../../raw_data/"

p = subprocess.Popen("ls " + path_to_raw_data, stdout = subprocess.PIPE, shell = True)
(output, err) = p.communicate()

p_status = p.wait() # Get error information incase there is one

# Input preparation and pre-cleaning
list_o_files = str(output).split("\\n") # Split resulting lines at the 

if list_o_files[0][1:2] == "b'": 
	list_o_files[0] = list_o_files[0].replace("b'","") # Removes the "b" signifying the sample as a bit type
		
if list_o_files[-1] == "'": 
	list_o_files = list_o_files[1:-1]

# Set up dictionary for sample read_group successes! 
sample_rg_dict = {}
arb_cutoff = 20

# Intake actual fastq files from symlinks and read through lines
for filenames in list_o_files: 

	lineNum = 0 
	symlink_filename = os.readlink(path_to_raw_data + filenames)

	shorter_path = "/".join(path_to_raw_data.split("/")[0:2]) + "/" # Removes "raw_data" to map symlink to file accurately

	readGroup = ""

	with gzip.open(shorter_path + symlink_filename) as tfile:

		success_count = 0 

		for lines in tfile: 
			lines = str(lines) # Convert from bit type to string

			if lineNum == 0: #Extract the sample name to ID into dictionary
				sample_name = lines[3:13]
				readGroup = "ID:" + sample_name + " PU:" + sample_name + " LB:" + sample_name + " PL:ILLUMINA SM:" + sample_name
				# print(readGroup)
				sample_rg_dict[ sample_name ] = [readGroup, 1] # Default to automatically entering read groups

			if lineNum > arb_cutoff: # Arbitrary cutoff where we stop reading in lines to process read group
				break

			if (lineNum % 4) == 0: 
				if len(lines.split(":")) >= 9: 
					success_count += 1

			lineNum += 1

		if success_count >= arb_cutoff: 
			sample_rg_dict[ sample_name ][1] = 0 # If the criteria is met for all sample reads, turn off the read groupings

		
#	
# 	gzipped_fastq_locale = symlink_path

# 	lineNum = 0

# 	BIG_BOOL_LIST = []

# 	for line in tfile: 

# 		with gzip.open(path_to_raw_data + gzipped_fastq_locale) as fastq_file: 

# 			for line in fastq_file: 
# 				if lineNum > 400: 
# 					break
				
# 				if (lineNum % 4) == 0: 
# 					nine_on_the_line_bool = (len(str(line).split(":")) >= 9) # Check if there are there 9 items in the list
# 					BIG_BOOL_LIST.append(nine_on_the_line_bool)

# 				lineNum += 1

# print(sum(BIG_BOOL_LIST))


# Open the file; find header lines and get 'em counted. Make sure they are less than

