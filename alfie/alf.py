import argparse

import seqio

"""
command line executable script for processing data with alfie

TODO - build the fasta single batch process, then adapt and repeat for the other 3 situations.

globals 
- build a dict of the outfile names where the keys are the kingdom numeric classes
- need the model stored in the program and loaded for the testing.

"""

def main():
	parser  = argparse.ArgumentParser(prog = "alfie",
		description = "Alfie\n"+\
		"a command line tool for kingdom-level classification and processing of DNA\n"+\
		"in fasta or fastq format")
	parser.add_argument("-f", "--file", type = str,  
		help = "The file of input sequences to classify."+\
		"Input can be either fasta or fastq formatfile type inferred from the extension.\n"+\
		"fasta: '.fasta' or '.fa'\n"+\
		"fastq: '.fastq' or '.fq'\n")
	parser.add_argument("-b", "--batch", type = int , default = 0, 
		help = "should the input file be processed in batches ofsequences?"+\
		"Default is False, passing an integer indicating the batch size to this flag"+\
		"will enable sub batches and decrease"+\
		"the amount of data stored in ram at one time. Tradeoff is slower processing")

	args = parser.parse_args()

	#check if fasta or fastq input
	ftype = seqio.file_type(args.file)

	# build the output filenames
	kingdom_outfiles = outfile_dict(args.file)

	if ftype == 'fasta': 
		if args.batch == 0:
			# batch fasta processing

			

		else:
			# full file fasta processing
			seq_records = seqio.read_fasta(args.file)

			# once read in generate the kmer data
			for entry in seq_records:
				entry['kmer_data'] = KmerFeatures(entry['name'], entry['sequence'])

			# turn kmer data into numpy array of proper structure
			# get the kmer_freqs for each 
			# kmer_data entry

			# load the model

			# pass the array of kmer data entries to the model

			# use the predictions to write the entries to the correct file
			# need the numeric encoding of the kingdoms stored here.



	elif ftype == 'fastq':
		if args.batch == 0:
			# batch fastq processing

		else:
			# full file fastq processing
