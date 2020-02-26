"""
note putting here in case I want to make a training api as well, if this is the case
then a sklearn dependency is introduced.
"""

import numpy as np

import alfie.seqio as seqio
from alfie.kmerseq import KmerFeatures


def upsample_fragments(seq_df, label_col = 'kingdom',  kmers=[4], **kwargs):
	"""
	take in a dataframe with dna in a column names 'clean_dna', conduct repeated
	subsampling of the sequence and generate kmer information for each subsample.

	- you can control the subsampling protocol by passing seq_samples keyword arguments
	to this function.

	returns a dict of lists: 
	keys(ids, label, data)

	ids - the processids
	label - the kingdom column by default, this can be changed as needed 
	data - the kmer array is 6mer+3mer as per previous
	"""

	#stores tuples of (processid, kingdom, kmer_freqs)
	samples = {'ids': [],
				'labels': [],
				'data': [],
				'seq': []}
	# data are split to train and test
	# now do the upsampling and generation of the output data files.
	for entry in seq_df.iterrows():
		processid =  entry[1].processid
		label = entry[1][label_col]
		seq =  entry[1].clean_dna

		sub_seqs = seq_samples(seq, **kwargs)

		for s in sub_seqs:
			k_seq = KmerFeatures(processid, s, kmers=kmers)
			samples['ids'].append(processid)
			samples['labels'].append(label)
			samples['data'].append(k_seq.kmer_freqs)
			samples['seq'].append(s)

	return samples


def shuffle_unison(x, y):
	"""
	shuffle the X and y numpy arrays in unison
	should be used if you're upsampling wiht the upsample_fragments function
	"""
	assert len(x) == len(y)
	p = np.random.permutation(len(x))
	return x[p], y[p]


def stratified_train_test_valid_taxon(input_data, class_col, test_size = 0.3):
	"""
	take a dataframe and conduct Stratified a train/test split based of a 
	user defined categorical column
	"""

	#split off a test/valid set, 30% of the data total
	strat_index = StratifiedShuffleSplit(n_splits=1, test_size=test_size, random_state=1738)

	for train_index, test_valid_index in strat_index.split(input_data, input_data[class_col]):
		X_train, X_test = input_data.loc[train_index], input_data.loc[test_valid_index] 


	return X_train, X_test

