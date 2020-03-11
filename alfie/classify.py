"""
A module with functions for sequence classification and decoding of classifications.

==========
Functions
==========

classify_records - Classify a series of DNA sequence records with the designated neural network.

decode_predictions - Decode numeric predictions to strings.

"""
import numpy as np

from alfie import dnn_k_four
from alfie.kmerseq import KmerFeatures


def classify_records(seq_records, dnn_model = dnn_k_four, k = 4):
	"""
	Classify a series of DNA sequence records with the designated neural network.

	The function takes a set of input sequences in the format returned by the seqio module's
	read_fasta or reads fastq functions. The `classify_records` function takes the 
	input sequences, generates a kmerseq.Kmerfeatures class instance for each
	sequence and passes the feature through the specified dnn_model to obtain a prediction. 

	Arguments
	---------
	seq_records : list, a list of sequence records. Where each record is a 
		dictionary with the keys 'name' (identifying string - header line) and 
		'sequence' (the sequence line of the fasta entry). Other keys permitted but unused

	dnn_model : tensorflow_model, A sequential tensorflow neural network. By default the interna
		kingdom-level classifier model is used. A user may specify a custom model, if the custom
		model utilizes a different kmer feature size, the k parameter must be altered accordingly. 

	k : int, the kmer input feature sizes corresponding to the dnn_model being uses. The
		features generated for each record will be kmer frequencies for size k.

	Returns
	---------

	out1, out2 : (list, array) out1 is a list of sequence records, with the new key, value
		pair 'kmer_data' added to each record's dictionary. out2 is an array classifications,
		whose length corresponds to the the length of the list of sequence records.

	Examples
	---------
	
	#load the example data from alfie
	>>> from alfie import example_fasta
	
	>>> seq_records, predictions = classify_records(example_fasta)
	
	>>> seq_records[0].keys()
	dict_keys(['name', 'sequence', 'kmer_data'])
	
	# position 0 of the predictions is the prediction for record 1 in the seq_records list
	# the prediction is a number, corresponding to the classes predicted by the given model
	>>> predictions[0]
	3
	"""


	for entry in seq_records:
		entry['kmer_data'] = KmerFeatures(entry['name'], entry['sequence'], k=k)

	vals = np.array([seq_records[i]['kmer_data'].kmer_freqs for i in range(len(seq_records))])
	
	yht_out = dnn_model.predict(vals)

	predictions = np.argmax(yht_out, axis = 1)

	return seq_records, predictions


def decode_predictions(predictions,
						tax_list = ["animalia","bacteria","fungi","plantae","protista",]):
	"""
	Decode numeric predictions to strings.

	Take in an array of numeric predictions from classify_records() and 

	Default decoding to kingdom. A custom tax_list with the classifications corresponding
	to a custom neural network's numeric precidions can be procided as well.
	
	Note: it is best practice for encoding and decoding to treat labels in alphabetical order.
	(this is the sklearn default)
	
	Arguments
	---------
	predictions : list like object, a list of numeric encoded predictions.

	tax_list : list, a list of strings indicating what the numeric predictions should
		be decoded to. By default kingdom labels are utilized.

	Returns
	---------
		out : list, a list of decoded strings.
	
	Examples
	---------
	#load example data
	>>> from alfie import example_fasta
	#generate predictions for examples, using default kingdom model
	>>> seq_records, predictions = classify_records(example_fasta)
	# decode predictions
	>>> predicted_kingdoms = decode_predictions(predictions)
	
	>>> predicted_kingdoms[:5]
	['plantae', 'bacteria', 'protista', 'animalia', 'animalia']
	"""
	class_dict = {}
	for i, x in enumerate(tax_list):
		class_dict[i] = x

	outpred = [class_dict[i] for i in predictions]

	return outpred

