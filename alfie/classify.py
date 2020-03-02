
import tensorflow
import numpy as np

from alfie import dnn_k_four
from alfie.kmerseq import KmerFeatures


def classify_records(seq_records, dnn_model = dnn_k_four, k = [4]):
	"""
	take a series of DNA sequence records and classify them with the 
	"""

	for entry in seq_records:
		entry['kmer_data'] = KmerFeatures(entry['name'], entry['sequence'], kmers=k)

	vals = np.array([seq_records[i]['kmer_data'].kmer_freqs for i in range(len(seq_records))])
	
	yht_out = dnn_model.predict(vals)

	predictions = np.argmax(yht_out, axis = 1)

	return seq_records, predictions