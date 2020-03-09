"""
Sequence input and output functions.

This module contains functions for the reading and writing of sequence data from fasta or fastq
formatted files. 

==========
Input functions
==========

read_fasta : 
read_fastq : 

iter_read_fasta :
iter_read_fastq :


==========
Output functions
==========

write_fasta :
write_fastq :

==========
Support functions
==========

file_type : Take in a filename and determine if the extension indicates a fasta or fastq file.
outfile_dict :
process_fastq_record :

"""

import os 
import copy

def file_type(s):
	"""
	Take in a filename and determine if the extension indicates a fasta or fastq file.

	Arguments
	---------
	s : str, a filename string

	Returns
	---------
	out : string, either 'fasta' or 'fastq' if file has an accepted extension, or ValueError.

	Examples
	---------
	>>> file_type("example_file.fasta")
	"fasta"
	>>> file_type("example_file.fa")
	"fasta"
	>>> file_type("example_file.fastq")
	"fastq"
	>>> file_type("example_file.fq")
	"fastq"
	>>> file_type("example_file.txt")
	ValueError: Input file must be in fasta or fastq format. Accepted file extensions: fa, fq, fasta, or fastq.	

	"""
	suffix = s.split(".")[-1]

	if suffix == "fa" or suffix == "fasta":
		return "fasta"

	elif suffix == "fq" or suffix == "fastq":
		return "fastq"

	else:
		raise ValueError("Input file must be in fasta or fastq format. "+\
			"Accepted file extensions: fa, fq, fasta, or fastq.")


def outfile_dict(filename, 
	labels = ["animalia", "bacteria", "fungi", "plantae", "protista"],
	folder_prefix = "alfie_out/"):
	""" 
	Build a dictionary of output filenames for classified sequences.

	Arguments
	---------
	filename - str, 
	labels - list, the list of labels (str) corresponding to the classification encoding. By default
		kingdom labels are utilized. To interface with sklearn's LabelBinarizer, labels should
		be in alphabetical order. 
	folder_prefix - str, the desired output location, prefix added to the output filenames. 
		By default, a new folder named 'alfie_out/' is generated. Passing 'folder_prefix = None'
		will omit the prefix, and files will be output to the current working directory and no
		new folder will be generated.

	Returns
	---------
	out : dict, keys are numbers denoting the order of the input labels and the values are
		a strings in the format: folder_prefix + label + filename

	Examples
	---------
	>>> outfile_dict('test_file.fasta', folder_prefix = None)
	{0: 'animalia_test_file.fasta',
	 1: 'bacteria_test_file.fasta',
	 2: 'fungi_test_file.fasta',
	 3: 'plantae_test_file.fasta',
	 4: 'protista_test_file.fasta'}

	>>> outfile_dict('test_file.fastq', labels = ['hot_dog','not_hot_dog'], folder_prefix = None)
	{0: 'hot_dog_test_file.fastq', 1: 'not_hot_dog_test_file.fastq'}

	"""
	f_stripped = filename.split('/')[-1]

	if folder_prefix != None:
		if os.path.isdir(folder_prefix) == False:
			os.mkdir(folder_prefix)
	else:
		folder_prefix = ''

	k_files = {}

	for i, x in enumerate(labels):
		k_files[i] = folder_prefix + x + "_" + f_stripped

	return k_files



def read_fasta(filename):
	""" 
	A read function that takes the data in from a fasta file.
	Outputs the data to a list of records in Bio.SeqRecord format


	Arguments
	---------

	filename

	Returns
	---------

	Examples
	---------
	# load the path to the alfie example file
	>>> from alfie import ex_fasta_file

	>>> data = read_fasta(ex_fasta_file)
	#data are a list of dictionaries with the keys 'name' and 'sequence'
	>>> data[0].keys()
	dict_keys(['name', 'sequence'])
	
	"""
	seq_records = []

	record = {"name" : None, "sequence" : ""}

	with open(filename) as file:
		for line in file:
			#if we hit a new record
			if line[0] == ">":
				#if current record, append to the record list
				if record["name"] != None:
					seq_records.append(copy.copy(record))	
				record["name"] = line[1:].rstrip()
				record["sequence"] = ""
			else:
				record["sequence"] += line.rstrip()

	seq_records.append(record)	

	return seq_records


def iter_read_fasta(filename, batch = 1000):
	"""	
	iteratively read in a fasta file, yielding batches of sequences. 
	batch size defined as argument, default is 1000 sequences


	Arguments
	---------

	Returns
	---------

	Examples
	---------


	"""
	seq_records = []

	record = {"name" : None, "sequence" : ""}

	with open(filename) as file:
		for line in file:
			#if we hit a new record
			if line[0] == ">":
				#if current record, append to the record list
				if record["name"] != None:
					seq_records.append(copy.copy(record))
					
					if len(seq_records)	== batch:
						yield seq_records
						seq_records = []

				record["name"] = line[1:].rstrip()
				record["sequence"] = ""
			else:
				record["sequence"] += line.rstrip()

	seq_records.append(record)	

	yield seq_records



def process_fastq_record(lines):
	""" 
	Take the four lines of a fastq record and create a dictonary for the record


	Arguments
	---------

	Returns
	---------

	Examples
	---------

	"""
	ks = ['name', 'sequence', 'strand', 'quality']
	return {k: v for k, v in zip(ks, lines)}


def read_fastq(filename):
	""" 
	This function takes a fastq filename and will read in the records in the file,
	constructing a dictonary for each record with the keys: 'name', 'sequence', 'optional', 'quality'.
		
	A list containing the dictonaries for all of the records will be returned.		

	Arguments
	---------

	Returns
	---------

	Examples
	---------

	"""
	records = []
	n = 4
	with open(filename, 'r') as file:
		lines = []
		for line in file:
			lines.append(line.rstrip())
			if len(lines) == n:
				record = process_fastq_record(lines)
				records.append(record)
				lines = []

	return records


def iter_read_fastq(filename, batch = 1000):
	""" 
	iteratively read in a fastq file, yielding batches of sequences. 
	batch size defined as argument, default is 1000 sequences

	Arguments
	---------

	Returns
	---------

	Examples
	---------

	"""
	records = []
	n = 4

	with open(filename, 'r') as file:
		lines = []
		for line in file:
			lines.append(line.rstrip())
			if len(lines) == n:
				record = process_fastq_record(lines)
				records.append(record)
				lines = []

				if len(records) == batch:
					yield records
					records = []
				
	if records:
		yield records


def write_fasta(entry, filename, append_seq = True):
	"""
	Takes a sequence record or list of sequence records. 
	Sequence record in dictonary format with the keys: 'name', 'sequence'
	Additional keys are permitted but unused.

	Will write the sequence to file, by default the sequence is appended to exiting entries,
	passing append_seq = False will overwrite the current data.



	Arguments
	---------

	Returns
	---------

	Examples
	---------

	"""
	if file_type(filename) != 'fasta':
		raise ValueError("output file does not have fasta extension.")

	if type(entry) == dict:
		entry = [entry]

	outstring = ''

	for x in entry:
		str_x = f">{x['name']}\n{x['sequence']}\n"
		outstring+=str_x

	if append_seq == True:
		mode = "a"
	else:
		mode = "w"
	file = open(filename, mode)
	file.write(outstring)
	file.close()


def write_fastq(entry, filename, append_seq = True):
	"""
	Takes a sequence record or list of sequence records. 
	Sequence record in dictonary format with the keys: 'name', 'sequence', 'plus', 'quality'
	Additional keys are permitted but unused.

	Will write the sequence to file, by default the sequence is appended to exiting entries,
	passing append_seq = False will overwrite the current data.


	Arguments
	---------

	Returns
	---------

	Examples
	---------

	"""
	if file_type(filename) != 'fastq':
		raise ValueError("output file does not have fastq extension.")

	if type(entry) == dict:
		entry = [entry]

	outstring = ''

	for x in entry:
		str_x = f"{x['name']}\n{x['sequence']}\n{x['strand']}\n{x['quality']}\n"
		outstring+=str_x

	if append_seq == True:
		mode = "a"
	else:
		mode = "w"

	file = open(filename, mode)
	file.write(outstring)
	file.close()


