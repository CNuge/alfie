"""
A module with functions to aid in training of custom alignment-free DNA classifiers.

==========
Functions
==========

alfie_dnn_default : Construct a simple neural network for sequence classification. 

process_sequences : Conduct subsampling of the sequences and generate kmer information for sequence.

sample_seq : Take a full sequence and return a list of random subsamples.

shuffle_unison : Shuffle the two input numpy arrays in unison.

stratified_taxon_split : Conduct a stratified train/test split based on a user defined categorical column.

"""

import numpy as np
import pandas as pd
import tensorflow as tf
from sklearn.model_selection import StratifiedShuffleSplit

from alfie.kmerseq import KmerFeatures


def stratified_taxon_split(input_data, class_col, test_size = 0.3, silent = False):
	"""
	Conduct a stratified train/test split based on a user defined categorical column.
	
	Arguments
	---------
	
	input_data : pandas.DataFrame, a dataframe to be split into a train and test set.

	class_col : string, the column of the input data with the categories to stratify 
		between the train and test set.

	test_size : double, the proportion of the input data to be included in the test split. 

	silent : bool, should the split criteria be echoded, defualt is True.

	Returns
	---------
	out1, out2 : pandas.DataFrame, out1 is the training data frame, out2 is the test data frame.

	Examples
	---------
	#initiate a similated dataframe
	>>> data = pd.DataFrame({"phylum" : ["Mollusca"]*10 + ["Arthropoda"] * 15,
	>>>						"data_col" : [np.random.randint(100) for x in range(25)]})
	#split on the column phylum, contians the classifications
	>>> train, test = stratified_taxon_split(data, class_col = "phylum", 
	>>>			test_size = .2, silent = True)
	# 80% of data in train
	>>> train.shape
	# index order is randomized
	>>> train.index
	# 20% of data in test
	>>> test.shape
	"""
	if silent == False:
		print(f'Conducting train/test split, split evenly by: {class_col}')

	#split off a test/valid set, 30% of the data total
	strat_index = StratifiedShuffleSplit(n_splits=1, test_size=test_size, random_state=1738)

	for train_index, test_valid_index in strat_index.split(input_data, input_data[class_col]):
		train, test = input_data.loc[train_index], input_data.loc[test_valid_index] 


	return train, test


def sample_seq(seq, min_size = 200, max_size = 600, n = 1, seed = None):
	"""
	Take a full sequence or list of sequences and return a list of random subsamples.

	Samples will be of a random length subset of the input seq. The min and max size of
	the random subset are defined by the min_size and max_size parameters. 

	Arguments
	---------
	seq : string or list, the sequence, or list of sequences, to randomly subsample.

	min_size : int, the minimum size of the random subsample. Default is 200.

	max_size : int, the maximum size of the random subsample. Default is 600.

	n : int, the number of random samples to generte from each input sequence.
		Default is 1 (no upsampling).

	seed : int, a random seed for repeatable random sampling.

	Returns
	---------

	Examples
	---------


	"""
	#list of output sequences
	outseqs = []
	#set random seed if passed by user
	if seed != None:
		random.seed(seed)
	#set the max to seq length if its shorter
	if max_size > len(seq):
		max_size = len(seq)
	if min_size > len(seq):
		raise ValueError("Minimum sample size exceeds sequence length")
	#get the set of random window sizes
	win_sizes = [np.random.randint(min_size, max_size) for x in range(n)]
	#for each window size, randomly subset the sequence by choosing a start point
	#and slicing the seq.
	for win_x in win_sizes:
		win_start = np.random.randint(0, (len(seq) - win_x))
		subseq = seq[win_start:(win_start+win_x)]
		outseqs.append(subseq)

	return outseqs


def process_sequences(seq_df, id_col = 'processid',
							seq_col = 'sequence', 
							label_col = 'kingdom',  k=4, **kwargs):
	"""
	Conduct subsampling of the sequences and generate kmer information for sequence.


	take in a dataframe with dna sequence, label and id information.

	returns a dict of lists: 
	keys(ids, label, data, seq)

	ids - the sequence IDs
	label - the sequence label column 
	data - the kmer array frequency as per previous
	seq - the DNA sequence used to generate the kmer array
	
	Arguments
	---------

	Returns
	---------

	Examples
	---------


	"""

	#stores tuples of (processid, kingdom, kmer_freqs)
	samples = {'ids': [],
				'labels': [],
				'data': [],
				'seq': []}
	# data are split to train and test
	# now do the upsampling and generation of the output data files.
	for entry in seq_df.iterrows():
		processid =  entry[1][id_col]
		label = entry[1][label_col]
		seq =  entry[1][seq_col]

		sub_seqs = sample_seq(seq, **kwargs)

		for s in sub_seqs:
			k_seq = KmerFeatures(processid, s, k=kmers)
			samples['ids'].append(processid)
			samples['labels'].append(label)
			samples['data'].append(k_seq.kmer_freqs)
			samples['seq'].append(s)

	return samples


def shuffle_unison(x, y):
	"""
	Shuffle the two input numpy arrays in unison.

	Should be used if you're upsampling with the upsample_fragments function
	

	Arguments
	---------

	Returns
	---------

	Examples
	---------

	"""
	assert len(x) == len(y)
	p = np.random.permutation(len(x))
	return x[p], y[p]


def alfie_dnn_default(hidden_sizes = [100], dropout = 0.2,
						in_shape = 256, n_classes = 5):
	"""
	Construct a simple neural network for sequence classification. 


	Arguments
	---------
	hidden_sizes - neuron sizes for the hidden layers
				n_hidden is implict param - equal to the length of hidden layers list
	dropout - dropout applied after each hidden layer, for no dropout pass 0 
	in_shape - the number of predictors this is for 1d inputs
	n_classes - the number of output classes

	Returns
	---------
	out : 

	Examples
	---------


	"""
	#initiate the model
	model = tf.keras.models.Sequential()
	#specify the in layer, denoting size
	model.add(tf.keras.layers.Dense(100, input_shape=(in_shape,) , activation = 'relu'))

	n_hidden = len(hidden_sizes)
	
	for i in range(0,n_hidden):
		model.add(tf.keras.layers.Dense(hidden_sizes[i], activation = 'relu'))
		if dropout != 0:
			model.add(tf.keras.layers.Dropout(dropout))


	model.add(tf.keras.layers.Dense(n_classes, activation = 'softmax'))

	model.compile(loss = 'sparse_categorical_crossentropy', 
						optimizer = 'adam', 
						metrics = ['accuracy'] )
	
	return model
