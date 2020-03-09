"""
alfie: Alignment-free identification of eDNA

==========
Data and models
==========

The alfie package contains the following data structures that can be imported with the syntax: 
from alfie import ___

ex_fasta_file : str, a path to a fastq file with 100 example COI-5P DNA barcode sequences

ex_fastq_file : str, a path to a fasta file with 100 example COI-5P DNA barcode sequences

example_fasta : list, a list of DNA sequences records. Each record contains data 
	derieved from a fasta file in dictionary format.

example_fastq : list, a list of DNA sequences records. Each record contains data 
	derieved from a fastq file in dictionary format.

dnn_k_four : a pre-trained tensorflow neural network for alignment-free, kingdom-level classification
	of COI-5P barcode sequences. The model is designed to take 4-mer frequencies as input and make a 
	multiclass prediction of the kingdom of origin for the input. Both the input 4-mers and the output
	kingdom-level classification are encoded values in alphabetical order. This is the defualt neural
	network deployed by the alfie command line interface, and by the 'classify_records' function in the
	alfie module 'classify'.

dnn_k_six : a pre-trained tensorflow neural network for alignment-free, kingdom-level classification
	of COI-5P barcode sequences. The model is designed to take 6-mer frequencies as input and make a 
	multiclass prediction of the kingdom of origin for the input. Both the input 6-mers and the output
	kingdom-level classification are encoded values in alphabetical order.
"""

import os

import alfie.seqio as seqio
from tensorflow.keras.models import load_model

location = os.path.dirname(os.path.realpath(__file__))

ex_fasta_file = os.path.join(location, 'data', 'example_data.fasta')
ex_fastq_file = os.path.join(location, 'data', 'example_data.fastq')

fourmer_model_file = os.path.join(location, 'data', 'dnn_model_4mers')
sixmer_model_file = os.path.join(location, 'data', 'dnn_model_6mers')

example_fasta = seqio.read_fasta(ex_fasta_file)
example_fastq = seqio.read_fastq(ex_fastq_file)

dnn_k_four = load_model(fourmer_model_file)
dnn_k_six = load_model(sixmer_model_file)