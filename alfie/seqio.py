

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
			"Accepted file extensions: fa, f1, fasta, or fastq ")