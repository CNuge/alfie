import pytest
from alfie import classify

from alfie import example_fasta
	
def test_classification_worklow():

	#go through the fasta, extract actual classes stored in header
	predictions_expected = []

	for x in example_fasta:
		name = x["name"]
		str_class = name.split('_')[-1]

		predictions_expected.append(str_class)

	#make predictions on the decod
	seq_records, predictions = classify.classify_records(example_fasta)
	
	assert list(seq_records[0].keys()) == ["name", "sequence", "kmer_data"]

	predictions_actual = classify.decode_predictions(predictions)
	
	for i, x in enumerate(predictions_expected):
		assert x == predictions_expected[i]




def test_decode_predictions():

	tax_classes = ["bill", "george", "sue",]

	num_input = [0 , 1 , 2 , 2 , 1, 0,]

	expected_output = ["bill", "george", "sue", "sue", "george", "bill",]

	predictions_custom = classify.decode_predictions(num_input, 
														tax_list = tax_classes)

	for i, x in enumerate(predictions_custom):
		assert x == expected_output[i]
