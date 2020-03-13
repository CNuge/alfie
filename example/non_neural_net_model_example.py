"""
alfie can be used with non-neural network models as well.

Here this is demonstrated by training a support vector machine 
on the same set of data used in the example jupyter notebook (part 2),

"""

import numpy as np
import pandas as pd

from alfie.classify import classify_records, decode_predictions
from alfie.training import stratified_taxon_split, process_sequences

from sklearn.svm import LinearSVC
from sklearn.metrics import accuracy_score
from sklearn.preprocessing import LabelBinarizer


#####
# load the example data
#####
data = pd.read_csv('alfie_small_train_example.tsv', sep = '\t')


#####
# conduct the train test split
#####
train, test = stratified_taxon_split(data, class_col = 'class', test_size = 0.3, )

train_kmer_data = process_sequences(train, label_col = 'class')
test_kmer_data = process_sequences(test, label_col = 'class')

#####
# encode the predictor and response data
#####
print("building X arrays")
X_train = np.array(train_kmer_data['data'])
X_test = np.array(test_kmer_data['data'])

print("building y arrays")
y_train_raw =  train_kmer_data['labels']
y_test_raw = test_kmer_data['labels']

tax_encoder = LabelBinarizer()

y_train = tax_encoder.fit_transform(y_train_raw)
y_test = tax_encoder.transform(y_test_raw)


y_train = np.reshape(y_train, len(y_train))
y_test = np.reshape(y_test, len(y_test))


#####
# train a demo SVM model
#####

svm_params = {'C': 100.0, 'loss': 'squared_hinge', 'max_iter' : 10000}

svm_ann_demo = LinearSVC(**svm_params) 
svm_ann_demo.fit(X_train, y_train) 


#####
# generate simulated fasta input from the test data
#####

#new list of dictionaries
test_simulated_fasta = []

#only doing first 10 rows
for i in range(0, len(y_test)):
    x = test.iloc[i]
    #make new dictionary entry for the sequence
    new_record = {'name' : x['processid'], 
                 'sequence': x['sequence']}
    #append to the list
    test_simulated_fasta.append(new_record)

#observe the data format
test_simulated_fasta[0]


#####
# use the model to make predictions via the classify records function
#####

# note the parameter argmax = False, this is passed because the LinearSVC already
# outputs non-one hot encoded outputs.
test_out, test_predictions = classify_records(test_simulated_fasta, 
												model = svm_ann_demo, argmax = False)

#evaluate against the true values
dnn_score = accuracy_score(y_test, test_predictions)

print("accuracy of model on classifying the test data:")
print(dnn_score)

