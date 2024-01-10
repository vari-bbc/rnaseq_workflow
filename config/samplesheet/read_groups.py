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
	list_o_files = list_o_files[1:-1] # Take everything except the last thing.


# Set up dictionary for sample read_group successes! 
sample_rg_dict = {}
arb_cutoff = 20

symlinks_present_bool = 0

# Intake actual fastq files from symlinks and read through lines
for filenames in list_o_files: 

	lineNum = 0 
	aug_path = "/".join(path_to_raw_data.split("/"))  # Removes "raw_data" to map symlink to file accurately
	sample_name	= filenames.split("_")[0]
	# sample_name = "_".join(filenames.split(".")[0].split("_")[:-1]) # Full file name up to "_R1" or "_R2"

	try:
		symlink_filename = os.readlink(path_to_raw_data + filenames)
		final_path = aug_path + symlink_filename
		symlinks_present_bool = 1

	except FileNotFoundError: 
		print("Cannot find FastQ files or symlinks in the", path_to_raw_data, "directory. Make sure they are present.")
		quit()

	except OSError: 
		if symlinks_present_bool: # This will only be triggered if files successfully went through the os.readlink
			print("\nWarning: FastQ files may be mixed with symlinks of FastQ files in '" + path_to_raw_data + "'. Careful!\n" )
			symlinks_present_bool = 0 # Only prints once!

		final_path = path_to_raw_data + filenames

	readGroup = ""

	with gzip.open(final_path) as tfile:

		success_count = 0 

		for lines in tfile: 
			lines = str(lines) # Convert from bit type to string

			if lineNum == 0: #Extract the sample name to ID into dictionary
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

	tfile.close()


# Next thing to do is to write the 
new_line_list = []

# Next bit might be better as a function. That way, you can check, modify, and pass path names without the rest of the loop
with open("units_template.tsv", "r") as ufile:

	header_bool = 1

	for lines in ufile: 

		if header_bool: 
			header = lines.strip() + "\n"
			header_bool = 0
			new_line_list.append(header)
			continue

		line_list = lines.split("\t")

		sample = line_list[0]
		group = line_list[1]
		fq1 = line_list[2]
		fq2 = line_list[3]

		# If there's something funky going on with read groups, it's probably here: 
		if sample_rg_dict[sample_name][1]:
			RG = sample_rg_dict[sample_name][0]
		else: 
			RG = ""

		new_line = "\t".join([sample, group, fq1, fq2, RG])
		new_line_list.append(new_line)

ufile.close()

# Rewrite units_template.tsv
new_file = open("units_template.tsv", "w")
# new_file = open("units.tsv", "w")

for lines in new_line_list: 

	print(lines.strip(), file = new_file)
	
new_file.close()




