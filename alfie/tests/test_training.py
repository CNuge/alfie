"""Unit tests for the module: alfie.training """

import pytest
import alfie.training as training
import numpy as np
import pandas as pd

#NOTE : I'm trying this in pytest as opposed to the unittest module, will see how it goes.

def test_split():
	"""Tests for the stratified_taxon_split function."""
	data = pd.DataFrame({"phylum" : ["Mollusca"]*10 + ["Arthropoda"] * 15,
							"data_col" : [np.random.randint(100) for x in range(25)]})
	#split on the column phylum, contians the classifications
	train, test = training.stratified_taxon_split(data, class_col = "phylum", 
				test_size = .2, silent = True, seed = 1738)
	# 80% of data in train
	assert train.shape == (20, 2)
	# index order is randomized
	assert list(train.index) == [16, 13, 0, 17, 5, 3, 10, 
									9, 18, 24, 23, 14, 2, 1, 
									20, 12, 19, 6, 4, 22]

	assert test.shape == (5, 2)
	assert list(test.index) == [15, 21, 7, 11, 8]


def test_sample_seq():

	in_seq = "AAAAAAAAAATTTTTTTTTTGGGGGGGGGGCCCCCCCCCCAAAAAAAAAATTTTTTTTTTGGGGGGGGGG"

	out1 = training.sample_seq(in_seq, min_size = 25, max_size = 70, seed = 1738)
	expected1 = ['GGGCCCCCCCCCCAAAAAAAAAATTTTTTT']
		
	assert out1 == expected1

	out2 = training.sample_seq(in_seq, min_size = 25, max_size = 70, n = 2, seed = 1738)
	expected2 = ['ATTTTTTTTTTGGGGGGGGGGCCCCCCCCC',
 					'TTTTTTTTGGGGGGGGGGCCCCCCCCCCAAAAAAAAAATTTTTTTTTTGGGG']

	assert out2 == expected2


def test_process_sequences():

	ex_dat = pd.DataFrame({"processid" : ["ex1", "ex2", "ex3", "ex4", "ex5",],
							"sequence" : ["AAAAAG" * 50 , "AAATAA" * 50, "AAGAAA" * 50, "TTTTAT" * 50, "TCTTCT" * 50],
							"kingdom" : ["animalia", "bacteria", "fungi", "plantae", "protista"]})

	#process the example data with defaults
	out_dat = training.process_sequences(ex_dat)

	#dict with 4 equal lenght lists
	assert list(out_dat.keys()) == ['ids', 'labels', 'data', 'seq']
	assert len(out_dat['ids']) == len(ex_dat['processid'])

	#different size k, turn off the subsampling, output a dataframe
	out_dat2 = training.process_sequences(ex_dat, k = 2, 
									to_dataframe = True, 
										subsample = False) 

	#query dataframe properties
	assert list(out_dat2.columns) ==['ids', 'labels', 'data', 'seq']
	assert np.all(out_dat2.ids == ex_dat.processid)
	assert out_dat2['data'][0].shape == (16,)



def test_shuffle_unison():

	x = np.array([[1,2],
					[3,4],
					[5,6],
					[7,8]])
	y = np.array([[1,2],
					[3,4],
					[5,6],
					[7,8]])

	new_x, new_y = training.shuffle_unison(x, y, seed = 1738)

	#is x the same as before shuffle_unison?
	assert np.all(new_x == x) == False
	#have x and y been shuffled in unison?
	assert np.all(new_x == new_y)

	with pytest.raises(ValueError):
		training.shuffle_unison(x, np.array([[1,1],[1,2]]), seed = 1738)

def test_alfie_dnn_default():

	model1 = training.alfie_dnn_default(hidden_sizes = [10,4], in_shape = 4, n_classes = 2)
	
	assert list(model1.input.shape) == [None, 4]
	
	assert list(model1.output.shape) == [None, 2]

	assert model1.trainable

