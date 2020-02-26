
import argparse
import numpy as np
import tensorflow as tf

from alfie import dnn_k_four

import alfie.seqio as seqio
from alfie.kmerseq import KmerFeatures

def process_records(seq_records, dnn_model = dnn_k_four , k = [4]):

	for entry in seq_records:
		entry['kmer_data'] = KmerFeatures(entry['name'], entry['sequence'], kmers=k)

	vals = np.array([seq_records[i]['kmer_data'].kmer_freqs for i in range(len(seq_records))])
	
	yht_out = dnn_model.predict(vals)

	predictions = np.argmax(yht_out, axis = 1)

	return seq_records, predictions

def main():
	parser  = argparse.ArgumentParser(prog = "alfie",
		description = "alfie:\n"+\
		"a command line tool for kingdom-level classification and processing of DNA\n"+\
		"in fasta or fastq format")
	parser.add_argument("-f", "--file", type = str,  
		help = "The file of input sequences to classify.\n"+\
		"Input can be either fasta or fastq formatfile type inferred from the extension.\n"+\
		"fasta: '.fasta' or '.fa' \n"+\
		"fastq: '.fastq' or '.fq' \n")
	parser.add_argument("-m", "--model", type = str, default = '4mer',
		help = "A file with a trained neural network to evaluate sequences." +\
		"If no model is specified, the default 4mer model is used"+\
		"Testing has shown the default 4mer model to be ~99.5 percent accurate."+\
		"Ensure the kmer size corresponds to the network input architeture" +\
		"other pre trained COI-5P kmer models available at:"+\
		"https://github.com/CNuge/alfie/models")
	parser.add_argument("-b", "--batch", type = int , default = 0, 
		help = "should the input file be processed in batches ofsequences?"+\
		"Default is False, passing an integer indicating the batch size to this flag"+\
		"will enable sub batches and decrease"+\
		"the amount of data stored in ram at one time. Tradeoff is slower processing")
	parser.add_argument("-k", "--kmer", type = int , default = 4, 
		help = "The kmer size used to evaluate sequences. Options 4mer (default) or 6mer "+\
		"The kmer features generated will correspond to the given size "+\
		"Testing has shown a 4mer model to be optimal")

	args = parser.parse_args()

	file = args.file
	model_file = args.model
	kmer = [args.kmer]
	batch = args.batch

	if file == None:
		raise ValueError("must specify an input data file with the flag -f")
	if model_file == None:
		raise ValueError("must specify an input model file with the flag -m")
	
	if model_file == '4mer':
		dnn_model = dnn_k_four
	else:
		# load the tensorflow model
		dnn_model = tf.keras.models.load_model(model_file)

	#model = '4mer'
	#file = '../data/example_data.fasta'
	#file = '../data/example_data.fastq'


	#check if fasta or fastq input
	ftype = seqio.file_type(file)

	# build the output filenames
	kingdom_outfiles = seqio.outfile_dict(file)


	if ftype == 'fasta': 
		if batch != 0:
			# batch fasta processing
			for b in seqio.iter_read_fasta(file, batch):

				seq_records, predictions = process_records(b, dnn_model, kmer)

				for i, entry in enumerate(seq_records):
					outfile = kingdom_outfiles[predictions[i]]
					seqio.write_fasta(entry, outfile)

		else:
			# full file fasta processing
			seq_records = seqio.read_fasta(file)
			
			seq_records, predictions = process_records(seq_records, dnn_model, kmer)

			for i, entry in enumerate(seq_records):
				outfile = kingdom_outfiles[predictions[i]]
				seqio.write_fasta(entry, outfile)


	elif ftype == 'fastq':
		if batch != 0:
			# batch fastq processing
			for b in seqio.iter_read_fastq(file, batch):

				seq_records, predictions = process_records(b, dnn_model, kmer)

				for i, entry in enumerate(seq_records):
					outfile = kingdom_outfiles[predictions[i]]
					seqio.write_fastq(entry, outfile)

		else:
			# full file fastq processing
			seq_records = seqio.read_fastq(file)

			seq_records, predictions = process_records(seq_records, dnn_model, kmer)

			for i, entry in enumerate(seq_records):
				outfile = kingdom_outfiles[predictions[i]]
				seqio.write_fastq(entry, outfile)


if __name__ == '__main__':
	main()