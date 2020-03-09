"""
Sequence input and output functions.

This module contains functions for the reading and writing of sequence data from fasta or fastq
formatted files. 

==========
Input functions
==========

read_fasta : Read data from a fasta file.

read_fastq : Read data from a fastq file.

iter_read_fasta : Iteratively read data from fasta file. 

iter_read_fastq : Iteratively read data from fastq file. 

==========
Output functions
==========

write_fasta : Write a sequence record, or list of records, to a file in fasta format.

write_fastq : Write a sequence record, or list of records, to a file in fastq format.

==========
Support functions
==========

file_type : Take in a filename and determine if the extension indicates a fasta or fastq file.

outfile_dict : Build a dictionary of output filenames for classified sequences.

process_fastq_record : Create a dictionary from a list of the four lines of a fastq record.

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
	Read data from a fasta file.
	
	Arguments
	---------
	filename : str, the path to a file in fasta format.

	Returns
	---------	
	out : list, a list of sequence records, where each record is a dictionary
	with the keys 'name' (identifying string - header line) and 
	'sequence' (the sequence line of the fasta entry).

	Examples
	---------
	# load the path to the alfie example file
	>>> from alfie import ex_fasta_file
	# read in the data
	>>> data = read_fasta(ex_fasta_file)
	# data are a list of dictionaries with the keys 'name' and 'sequence'
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
	Iteratively read data from fasta file. 
	
	Arguments
	---------
	filename : str, the path to a file in fasta format.

	batch : int, the number of sequence records to be returned in each batch.
		The default is 1000.

	Returns
	---------	
	out : generator, will yield lists of sequence records for the
	given batch size. 
	The lists contain sequence records in dictionary format,
	with the keys 'name' (identifying string - header line) and 
	'sequence' (the sequence line of the fasta entry).

	Examples
	---------
	# load the path to the alfie example file
	>>> from alfie import ex_fasta_file
	# read in the data
	>>> data = iter_read_fasta(ex_fasta_file, batch = 10)
	# data is a generator, get the next batch
	>>> x = next(data)
	# each record in the list is a dictionary
	>>> x[0]
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
	""" Create a dictionary from a list of the four lines of a fastq record."""
	ks = ['name', 'sequence', 'strand', 'quality']
	return {k: v for k, v in zip(ks, lines)}


def read_fastq(filename):
	""" 
	Read data from a fastq file.
	
	Arguments
	---------
	filename : str, the path to a file in fastq format.

	Returns
	---------	
	out : list, a list of sequence records, where each record is a dictionary
	with the keys: 'name' (identifying string - header line) and 
	'sequence' (the sequence line of the fastq entry), 'strand' ('+'' or '-''), 
	and 'quality' (the PHRED quality string).

	Examples
	---------
	# load the path to the alfie example file
	>>> from alfie import ex_fastq_file
	# read in the file
	>>> data = read_fastq(ex_fastq_file)
	# data are a list of dictionaries with the keys:
	# 'name', 'sequence', 'strand', and 'quality'
	>>> data[0].keys()
	dict_keys(['name', 'sequence', 'strand', 'quality'])
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
	Iteratively read data from fastq file. 
	
	Arguments
	---------
	filename : str, the path to a file in fastq format.

	batch : int, the number of sequence records to be returned in each batch.
		The default is 1000.

	Returns
	---------	
	out : generator, will yield lists of sequence records for the
	given batch size. 
	The lists contain sequence records in dictionary format,
	with the keys: 'name' (identifying string - header line) and 
	'sequence' (the sequence line of the fastq entry), 'strand' ('+'' or '-''), 
	and 'quality' (the PHRED quality string).

	Examples
	---------
	# load the path to the alfie example file
	>>> from alfie import ex_fastq_file

	>>> data = iter_read_fasta(ex_fastq_file, batch = 10)
	#data is a generator, get the next batch
	>>> x = next(data)
	#each record in the list is a dictionary
	>>> x[0]
	dict_keys(['name', 'sequence', 'strand', 'quality'])
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
	Write a sequence record, or list of records, to a file in fasta format.

	Arguments
	---------
	entry : list or dict, a sequence record or list of sequence records. Each sequence record 
		in dictionary format with the keys: 'name' (becomes header line of fasta), 
		'sequence'(becomes the second line of the fasta).
	
	filename : str, the path to the output file. Must have extension '.fasta' or '.fa'.

	append_seq : bool, indicate if sequence should be appended to existing data in the file.
		Default is True. If False, the existing file is overwritten where appliciable.

	Returns
	---------
	out : No return, data written to file.

	Examples
	---------
	# single sequence output
	>>> seq_record = {"name" : "example1", "sequence" : "ATGCATGC"}
	
	>>> write_fasta(seq_record, "example_out.fasta")

	# multi sequence output
	>>> seq_records = [{"name" : "example1b", "sequence" : "ATGCATGC"},
	>>>					{"name" : "example2", "sequence" : "GGGGGGGAAAAA"} ]
	
	# append_seq = False, therefore overwrite the previous output
	>>> write_fasta(seq_records, "example_out.fasta", append_seq = False)
	"""
	if file_type(filename) != 'fasta':
		raise ValueError("Specified output file does not have fasta extension.")

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
	Write a sequence record, or list of records, to a file in fastq format.

	Arguments
	---------
	entry : list or dict, a sequence record or list of sequence records. Each sequence record 
		in dictionary format with the keys: 'name', 'sequence', 'strand', and 'quality'.
		The values are added to the corresponding line of the fastq record.
	
	filename : str, the path to the output file. Must have extension '.fastq' or '.fq'.

	append_seq : bool, indicate if sequence should be appended to existing data in the file
		Default is True. If False, the existing file is overwritten where appliciable.

	Returns
	---------
	out : No return, data written to file.

	Examples
	---------
	# single sequence output
	>>> seq_record = {"name" : "example1", "sequence" : "ATGCATGC", 
	>>>					"strand" : "+", "quality" : "~~~~~~~~"}
	
	>>> write_fastq(seq_record, "example_out.fastq")

	# multi sequence output
	>>> seq_records = [{"name" : "example1b", "sequence" : "ATGCATGC", 
	>>>					"strand" : "+", "quality" : "~~~~~~~~"},
	>>>					{"name" : "example2", "sequence" : "GGGGGGGAAAAA", 
	>>>					"strand" : "+", "quality" : "~~~9~~~+~*~~"}]
	
	# append_seq = False, therefore overwrite the previous output
	>>> write_fastq(seq_records, "example_out.fastq", append_seq = False)
	"""
	if file_type(filename) != 'fastq':
		raise ValueError("output file does not have fastq extension.")

	if type(entry) == dict:
		entry = [entry]

	outstring = ''

	for x in entry:
		str_x = f"@{x['name']}\n{x['sequence']}\n{x['strand']}\n{x['quality']}\n"
		outstring+=str_x

	if append_seq == True:
		mode = "a"
	else:
		mode = "w"

	file = open(filename, mode)
	file.write(outstring)
	file.close()


