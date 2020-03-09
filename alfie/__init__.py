
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