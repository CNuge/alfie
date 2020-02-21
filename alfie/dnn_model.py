def dnn_class_model(hidden_sizes = [256,128,64,32,16], dropout = 0.3,
						in_shape = 256, n_classes = 5):
	"""
	builds a simple deep neural network using the keras wrapper.
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


def init_dnn(weight_file = './'):
	initial_model = dnn_class_model()

