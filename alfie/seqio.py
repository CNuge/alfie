import os 
import copy

def file_type(s):
	"""
	Take in a filename passed to alfie via the command line.
	Determine if the file type is fasta, fq or error
	"""

	suffix = s.split(".")[-1]

	if suffix == "fa" or suffix == "fasta":
		return "fasta"

	elif suffix == "fq" or suffix == "fastq":
		return "fastq"

	else:
		raise ValueError("input file must be in fasta or fastq format."+\
			"Accepted file extensions: fa, fq, fasta, or fastq ")


def outfile_dict(filename):
	""" 
	build a dictonary with the output filenames.
	dict keys are the numeric encodings of the kingdom names
	"""
	f_stripped = filename.split('/')[-1]

	if os.path.isdir("alfie_out") == False:
		os.mkdir("alfie_out")

	kingdom_files = {
		0: "alfie_out/animalia_",
		1: "alfie_out/bacteria_",
		2: "alfie_out/fungi_",
		3: "alfie_out/plantae_",
		4: "alfie_out/protista_",
	}

	for k, v in kingdom_files.items():
		v += f_stripped
		kingdom_files[k] = v		

	return kingdom_files



def read_fasta(filename):
	""" 
	A read function that takes the data in from a fasta file.
	Outputs the data to a list of records in Bio.SeqRecord format
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
	"""
	ks = ['name', 'sequence', 'strand', 'quality']
	return {k: v for k, v in zip(ks, lines)}


def read_fastq(filename):
	""" 
	This function takes a fastq filename and will read in the records in the file,
	constructing a dictonary for each record with the keys: 'name', 'sequence', 'optional', 'quality'.
		
	A list containing the dictonaries for all of the records will be returned.		
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


