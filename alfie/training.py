"""
A module with functions to aid in training of custom alignment-free DNA classifiers.

==========
Functions
==========

alfie_dnn_default : Construct a neural network for alignment-free classification. 

process_sequences : Conduct subsampling of the sequences and generate kmer information for sequence.

sample_seq : Take a full sequence and return a list of random subsamples.

shuffle_unison : Shuffle the two input numpy arrays in unison.

stratified_taxon_split : Conduct a stratified train/test split based on a user defined categorical column.

"""

import numpy as np
import pandas as pd

from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout

from sklearn.model_selection import StratifiedShuffleSplit

from alfie.kmerseq import KmerFeatures


def stratified_taxon_split(input_data, class_col, test_size = 0.3, silent = False, seed = None):
	"""
	Conduct a stratified train/test split based on a user defined categorical column.
	
	The stratified nature of this split ensures that the frequencies of different
	classes in the input dataframe are maintained in both the train and test sets.
	Even sampling to reduce potential source bias.

	Arguments
	---------
	
	input_data : pandas.DataFrame, a dataframe to be split into a train and test set.

	class_col : string, the column of the input data with the categories to stratify 
		between the train and test set.

	test_size : double, the proportion of the input data to be included in the test split. 

	silent : bool, should the split criteria be echoded, defualt is True.

	seed : int, a random seed for repeatable random sampling. Default is None.


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
	strat_index = StratifiedShuffleSplit(n_splits=1, test_size=test_size, random_state=seed)

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

	seed : int, a random seed for repeatable random sampling. Default is None.

	Returns
	---------
	out : list, a list of the output sequences

	Examples
	---------
	#example string is 70bp in length
	>>> in_seq = "AAAAAAAAAATTTTTTTTTTGGGGGGGGGGCCCCCCCCCCAAAAAAAAAATTTTTTTTTTGGGGGGGGGG"
	
	#take a single sample of the input, note the min_bp must be less than len(in_seq)
	>>> sample_seq(in_seq, min_size = 25, max_size = 70, seed = 1738)
	['GGGCCCCCCCCCCAAAAAAAAAATTTTTTT']

	#upsample the input, n subsamples returned
	>>> sample_seq(in_seq, min_size = 25, max_size = 70, n = 2, seed = 1738)
	['ATTTTTTTTTTGGGGGGGGGGCCCCCCCCC',
	 'TTTTTTTTGGGGGGGGGGCCCCCCCCCCAAAAAAAAAATTTTTTTTTTGGGG']
	"""
	#list of output sequences
	outseqs = []
	#set random seed if passed by user
	if seed != None:
		np.random.seed(seed)
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
							label_col = 'kingdom',
							k = 4, 
							to_dataframe = False, 
							subsample = True, 
							**kwargs):
	"""
	Conduct subsampling of the sequences and generate kmer information for sequence.

	This function executes a sequence processing pipeline to generate inputputs for
	training the neural network. It takes a dataframe of inputs with columns containing:
	sequences, ids, and labels (additional columns ignored). Sequences are subsampled
	(option can be turned off) and Kmer count frequencies for the given value of k are
	generated. Output is by default a dictonary of lists
	
	Arguments
	---------
	seq_df : pd.DataFrame, a data frame with columns containing 
		dna sequences, labels (classifications), and id information.
	
	id_col : string, the name of the input dataframe column that contains the
		sequence identifiers. These become the 'name' arguments for the 
		KmerFeatures class instances. Default is 'processid'.

	seq_col : string, the name of the input dataframe column that contains the
		DNA sequences . These become the 'sequence' arguments for the 
		KmerFeatures class instances. Default is 'sequence'.

	label_col : string, the column used to generate the 'label'	

	to_dataframe : bool, logical indicating if the output should be returned as a
		pandas DataFrame. Default is False - returned as a dictionary of lists.
	
	subsample : bool, logical indicating if the input sequences should be subsampled
		with the sample_seq function. Default is true. If false, kmer frequencies are
		based on the unaltered input sequences and no upsampling is performed.

	**kwargs : additional keyword arguments to be passed to the sample_seq function.
		See: alfie.training.sample_seq for a list of arguments.

	Returns
	---------
	out : dict of lists of equal size. The keys are: ids, label, data, seq. Each index
		position is an individual sequence observation. The output can optionall be 
		provided as a pandas dataframe as well.
		key descriptions:
			ids - the sequence IDs
			label - the sequence label column 
			data - the kmer array frequencies for the given sequence
			seq - the subsample of the DNA sequence used to generate the kmer frequencies

	Examples
	---------
	#build a dataframe of artifical data
	>>> ex_dat = pd.DataFrame({
	>>>		"processid" : ["ex1", "ex2", "ex3", "ex4", "ex5"],
	>>>		"sequence" : ["AAAAAG"*50, "AAATAA"*50, "AAGAAA"*50, "TTTTAT"*50, "TCTTCT"*50],
	>>>		"kingdom" : ["animalia", "bacteria", "fungi", "plantae", "protista"]})

	#process the example data with defaults
	>>> out_dat = process_sequences(ex_dat)

	#dict with 4 equal lenght lists
	>>> out_dat.keys()
	dict_keys(['ids', 'labels', 'data', 'seq'])
	>>> len(out_dat['ids']) == len(ex_dat['processid'])

	#different size k, turn off the subsampling, output a dataframe
	>>> out_dat2 = process_sequences(ex_dat, k = 2, 
	>>>								to_dataframe = True, 
	>>>								subsample = False) 

	>>> out_dat2.columns
	Index(['ids', 'labels', 'data', 'seq'], dtype='object')
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

		if subsample == True:
			sub_seqs = sample_seq(seq, **kwargs)
		else:
			sub_seqs = [seq]

		for s in sub_seqs:
			k_seq = KmerFeatures(processid, s, k=k)
			samples['ids'].append(processid)
			samples['labels'].append(label)
			samples['data'].append(k_seq.kmer_freqs)
			samples['seq'].append(s)

	if to_dataframe == True:
		return pd.DataFrame(samples)
	
	return samples


def shuffle_unison(x, y, seed = None):
	"""
	Shuffle the two input numpy arrays in unison.

	Intended be used if you're upsampling with the upsample_fragments function
	to randomize the dataframe orders.

	Arguments
	---------
	x : np.array, the first array to shuffle

	y : np.array, the second array to shuffle

	seed : int, a random seed for repeatable random sampling. Default is None.

	Returns
	---------
	out1, out2 : np.arrays with shapes respectively matching the shapes of the
		x and y inputs

	Examples
	---------
	#two arrays, with equal values
	>>> x = np.array([[1,2],
	>>>				[3,4],
	>>>				[5,6],
	>>>				[7,8]])
	>>> y = np.array([[1,2],
	>>>				[3,4],
	>>>				[5,6],
	>>>				[7,8]])

	>>> new_x, new_y = shuffle_unison(x, y, seed = 1738)

	#is x the same as before shuffle_unison?
	>>> np.all(new_x == x)
	False
	#have x and y been shuffled in unison?
	>>> np.all(new_x == new_y)
	True	

	"""
	if len(x) != len(y):
		raise ValueError("The input arrays do not have equal lengths.")
	if seed != None:
		np.random.seed(seed)

	p = np.random.permutation(len(x))
	return x[p], y[p]


def alfie_dnn_default(hidden_sizes = [100], dropout = 0.3,
						in_shape = 256, n_classes = 5):
	"""
	Construct a neural network for alignment-free classification.

	Arguments
	---------
	hidden_sizes - neuron sizes for the hidden layers
				n_hidden is implict param - equal to the length of hidden layers list
	dropout : float, fraction of dropout applied after each hidden layer, for no dropout pass 0.
		Default is 0.3. 
	in_shape : int, the number of predictor variables, assumes 1d inputs. 
		Default is 256 (4mer size).
	n_classes - int, the number of output classes. Default is 5 (kingdoms).

	Returns
	---------
	out : a tensorflow sequential neural network.

	Examples
	---------
	# construct a simple model
	# two hidden layers (10 and 4 neurons respectively)
	# takes 4 input values (i.e. a 1mer model) 
	# makes binary predictions 
	>>> model1 = alfie_dnn_default(hidden_sizes = [10,4], in_shape = 4, n_classes = 2)
	#4 inputs
	>>> model1.input.shape
	TensorShape([None, 4])
	#binary output
	>>> model1.output.shape
	TensorShape([None, 2])
	>>> model1.trainable
	True
	"""
	#initiate the model
	model = Sequential()
	#specify the in layer, denoting size
	model.add(Dense(100, input_shape=(in_shape,) , activation = 'relu'))

	n_hidden = len(hidden_sizes)
	
	for i in range(0,n_hidden):
		model.add(Dense(hidden_sizes[i], activation = 'relu'))
		if dropout != 0:
			model.add(Dropout(dropout))


	model.add(Dense(n_classes, activation = 'softmax'))

	model.compile(loss = 'sparse_categorical_crossentropy', 
						optimizer = 'adam', 
						metrics = ['accuracy'] )
	
	return model

