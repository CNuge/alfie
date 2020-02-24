
import argparse
import numpy as np
import tensorflow as tf

import seqio
from kmerseq import KmerFeatures



def process_records(seq_records, dnn_model, model, k):

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
	parser.add_argument("-b", "--batch", type = int , default = 0, 
		help = "should the input file be processed in batches ofsequences?"+\
		"Default is False, passing an integer indicating the batch size to this flag"+\
		"will enable sub batches and decrease"+\
		"the amount of data stored in ram at one time. Tradeoff is slower processing")
	parser.add_argument("-m", "--model", type = str , default = "4mer", 
		help = "The neural network used to evaluate sequences. Options 4mer (default) or 6mer "+\
		"The models use different kmer feqture sizes to classify the input sequences"+\
		"Testing has shown the 4mer modelto be ~99.5 percent accurate. Test show 6mer model is more"+\
		"accurate (~99.8 percent accuracy), but this comes with the tradeoff of increased processing time.")

	args = parser.parse_args()

	file = args.file
	model = args.model
	batch = args.batch

	#model = '4mer'
	#file = '../data/example_data.fasta'
	#file = '../data/example_data.fastq'


	#check if fasta or fastq input
	ftype = seqio.file_type(file)

	# build the output filenames
	kingdom_outfiles = seqio.outfile_dict(file)

	# load the tensorflow model
	if model == '4mer':
		dnn_model = tf.keras.models.load_model('dnn_alfie/alf_dnn.h5')
		k = [4]
	elif model == '6mer':
		dnn_model = tf.keras.models.load_model('dnn_alfie/dnn_model_6mers.h5')
		k = [6]
	else:
		raise ValueError("valid model choices are '4mer' or '6mer'")

	if ftype == 'fasta': 
		if batch != 0:
			# batch fasta processing
			for b in seqio.iter_read_fasta(file, batch):

				seq_records, predictions = process_records(b, dnn_model, model, k)

				for i, entry in enumerate(seq_records):
					outfile = kingdom_outfiles[predictions[i]]
					seqio.write_fasta(entry, outfile)

		else:
			# full file fasta processing
			seq_records = seqio.read_fasta(file)
			
			seq_records, predictions = process_records(seq_records, dnn_model, model, k)

			for i, entry in enumerate(seq_records):
				outfile = kingdom_outfiles[predictions[i]]
				seqio.write_fasta(entry, outfile)


	elif ftype == 'fastq':
		if batch != 0:
			# batch fastq processing
			for b in seqio.iter_read_fastq(file, batch):

				seq_records, predictions = process_records(b, dnn_model, model, k)

				for i, entry in enumerate(seq_records):
					outfile = kingdom_outfiles[predictions[i]]
					seqio.write_fastq(entry, outfile)

		else:
			# full file fastq processing
			seq_records = seqio.read_fastq(file)

			seq_records, predictions = process_records(seq_records, dnn_model, model, k)

			for i, entry in enumerate(seq_records):
				outfile = kingdom_outfiles[predictions[i]]
				seqio.write_fastq(entry, outfile)


if __name__ == '__main__':
	main()