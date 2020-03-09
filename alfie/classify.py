
import tensorflow
import numpy as np

from alfie import dnn_k_four
from alfie.kmerseq import KmerFeatures


def classify_records(seq_records, dnn_model = dnn_k_four, k = 4):
	"""	classify a series of DNA sequence records with the designated neural network."""
	for entry in seq_records:
		entry['kmer_data'] = KmerFeatures(entry['name'], entry['sequence'], kmers=k)

	vals = np.array([seq_records[i]['kmer_data'].kmer_freqs for i in range(len(seq_records))])
	
	yht_out = dnn_model.predict(vals)

	predictions = np.argmax(yht_out, axis = 1)

	return seq_records, predictions


def decode_predictions(predictions,
						tax_list = ["animalia","bacteria","fungi","plantae","protista",]):
	"""
	Decode the numeric predictions to a list of strings.

	Take in an array of numeric predictions from classify_records() and 

	Default decoding to kingdom. A custom tax_list with the classifications corresponding
	to a custom neural network's numeric precidions can be procided as well.
	
	Note: it is best practice for encoding and decoding to treat labels in alphabetical order.
	(this is the sklearn default)
	"""
	class_dict = {}
	for i, x in enumerate(tax_list):
		class_dict[i] = x

	outpred = [class_dict[i] for i in predictions]

	return outpred

