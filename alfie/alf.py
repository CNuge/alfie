import argparse


"""
command line executable script for processing data with alfie

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



	args = parser.parse_args()


	#check if fasta or fq


	#