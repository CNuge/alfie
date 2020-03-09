"""
A module with functions to aid in training custom alignment-free DNA classifiers


==========
Functions
==========

sample_seq

stratified_taxon_split

process_sequences

shuffle_unison

alfie_dnn_default

"""
import tensorflow as tf

import numpy as np

from alfie.kmerseq import KmerFeatures

from sklearn.model_selection import StratifiedShuffleSplit


def stratified_taxon_split(input_data, class_col, test_size = 0.3, silent = False):
	"""
	Conduct a stratified train/test split based of a user defined categorical column.
	"""
	if silent == False:
		print(f'Conducting train/test split, split evenly by: {class_col}')

	#split off a test/valid set, 30% of the data total
	strat_index = StratifiedShuffleSplit(n_splits=1, test_size=test_size, random_state=1738)

	for train_index, test_valid_index in strat_index.split(input_data, input_data[class_col]):
		X_train, X_test = input_data.loc[train_index], input_data.loc[test_valid_index] 


	return X_train, X_test


def sample_seq(seq, min_size = 200, max_size = 600, n = 1, seed = None):
	"""
	Take a full sequence and return a list of random subsamples.

	Samples will be of a random length form within the defined sizes of the 
	"""
	#list of output sequences
	outseqs = []
	#set random seed if passed by user
	if seed != None:
		random.seed(seed)
	#set the max to seq length if its shorter
	if max_size > len(seq):
		max_size = len(seq)
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
							label_col = 'kingdom',  kmers=4, **kwargs):
	"""
	Conduct subsampling of the sequences and generate kmer information for sequence.


	take in a dataframe with dna sequence, label and id information.

	returns a dict of lists: 
	keys(ids, label, data, seq)

	ids - the sequence IDs
	label - the sequence label column 
	data - the kmer array frequency as per previous
	seq - the DNA sequence used to generate the kmer array
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
			k_seq = KmerFeatures(processid, s, kmers=kmers)
			samples['ids'].append(processid)
			samples['labels'].append(label)
			samples['data'].append(k_seq.kmer_freqs)
			samples['seq'].append(s)

	return samples


def shuffle_unison(x, y):
	"""
	Shuffle the two input numpy arrays in unison.

	Should be used if you're upsampling wiht the upsample_fragments function
	"""
	assert len(x) == len(y)
	p = np.random.permutation(len(x))
	return x[p], y[p]


def alfie_dnn_default(hidden_sizes = [100], dropout = 0.2,
						in_shape = 256, n_classes = 5):
	"""
	Construct a simple neural network for sequence classification. 


	hidden_sizes - neuron sizes for the hidden layers
				n_hidden is implict param - equal to the length of hidden layers list
	dropout - dropout applied after each hidden layer, for no dropout pass 0 
	in_shape - the number of predictors this is for 1d inputs
	n_classes - the number of output classes
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
