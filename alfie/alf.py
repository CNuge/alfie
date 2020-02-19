import argparse

import seqio

"""
command line executable script for processing data with alfie

TODO - build the fasta single batch process, then adapt and repeat for the other 3 situations.

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
	parser.add_argument("-b", "--batch", action = "store_true", 
		help = "should the input file be processed in batches of 1000 sequences?"+\
		"Default is False, passing this flag will enable sub batches and decrease"+\
		"the amount of data stored in ram at one time. tradeoff is slower ")

	args = parser.parse_args()

	#check if fasta or fastq input
	ftype = seqio.file_type(args.file)

	if ftype == 'fasta': 
		if args.batch == True:
			# batch fasta processing


		else:
			# full file fasta processing

			# read in the file, store entries as tuples 

			# once read in generate the kmer data

			# turn kmer data into numpy array of proper structure (bind 3mer and 6mer, 
			# and bind rows to make shape) 

			# load the model

			# pass the array of kmer data entries to the model

			# use the predictions to write the entries to the correct file
			# need the numeric encoding of the kingdoms stored here.



	elif ftype == 'fastq':
		if args.batch == True:
			# batch fastq processing

		else:
			# full file fastq processing
