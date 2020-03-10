import unittest

import training
import pandas as pd

class TrainingTests(unittest.TestCase):

	def test_split(self):

		data = pd.DataFrame({"phylum" : ["Mollusca"]*10 + ["Arthropoda"] * 15,
								"data_col" : [np.random.randint(100) for x in range(25)]})
		#split on the column phylum, contians the classifications
		train, test = stratified_taxon_split(data, class_col = "phylum", 
					test_size = .2, silent = True)
		# 80% of data in train
		train.shape
		# index order is randomized
		train.index
		test.shape


	def test_sample_sequences(self):
		in_seq = "AAAAAAAAAATTTTTTTTTTGGGGGGGGGGCCCCCCCCCCAAAAAAAAAATTTTTTTTTTGGGGGGGGGG"
	
		sample_seq(in_seq, min_size = 25, max_size = 70, seed = 1738)
		['GGGCCCCCCCCCCAAAAAAAAAATTTTTTT']

		sample_seq(in_seq, min_size = 25, max_size = 70, n = 2, seed = 1738)
		['ATTTTTTTTTTGGGGGGGGGGCCCCCCCCC',
	 	'TTTTTTTTGGGGGGGGGGCCCCCCCCCCAAAAAAAAAATTTTTTTTTTGGGG']


	def test_process_sequences(self):

		ex_dat = pd.DataFrame({"processid" : ["ex1", "ex2", "ex3", "ex4", "ex5",],
								"sequence" : ["AAAAAG" * 50 , "AAATAA" * 50, "AAGAAA" * 50, "TTTTAT" * 50, "TCTTCT" * 50],
								"kingdom" : ["animalia", "bacteria", "fungi", "plantae", "protista"]})

		#process the example data with defaults
		out_dat = process_sequences(ex_dat)

		#dict with 4 equal lenght lists
		out_dat.keys()
		dict_keys(['ids', 'labels', 'data', 'seq'])
		len(out_dat['ids']) == len(ex_dat['processid'])

		#different size k, turn off the subsampling, output a dataframe
		out_dat2 = process_sequences(ex_dat, k = 2, 
										to_dataframe = True, 
											subsample = False) 

		out_dat2.columns
		Index(['ids', 'labels', 'data', 'seq'], dtype='object')

	def test_shuffle_unison(self):

		x = np.array([[1,2],
						[3,4],
						[5,6],
						[7,8]])
		y = np.array([[1,2],
						[3,4],
						[5,6],
						[7,8]])

		new_x, new_y = shuffle_unison(x, y, seed = 1738)

		#is x the same as before shuffle_unison?
		np.all(new_x == x)
		False
		#have x and y been shuffled in unison?
		np.all(new_x == new_y)


	def test_nn_constriction(self):

		dnn_1mer = training.alfie_dnn_default()
		model1 = alfie_dnn_default(hidden_sizes = [10,4], in_shape = 4, n_classes = 2)
		
		model1.input.shape
		TensorShape([None, 4])
		
		model1.output.shape
		TensorShape([None, 2])
		
		model1.trainable

