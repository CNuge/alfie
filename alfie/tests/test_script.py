import os
import sys
import pytest

from alfie import alf
from alfie import ex_fasta_file, ex_fastq_file


def test_argparser():

	#test no arguments - gives all defaults
	out1 = alf.alfie_parser([])

	assert out1.file == None
	assert out1.model == '4mer'
	assert out1.kmer == 4
	assert out1.batch == 0
	assert out1.classes == 'kingdoms'

	out2 = alf.alfie_parser(['-f', 'example.fasta',
							'-k', '10',
							'-m', 'wacky_custom'])

	assert out2.file == 'example.fasta'
	assert out2.model == 'wacky_custom'
	assert out2.kmer == 10
	assert out2.batch == 0
	assert out2.classes == 'kingdoms'


def test_main_with_args():

	#test empy args processed properly
	sys.argv = ['alfie']
	with pytest.raises(ValueError):
		alf.main()

	#no f flag value processed properly
	sys.argv = ['alfie', "-f", ""]
	with pytest.raises(ValueError):
		alf.main()

	#process a fasta file
	sys.argv = ['alfie', "-f", ex_fasta_file]
	alf.main()

	fasta_main_outputs = ['animalia_example_data.fasta',
							'bacteria_example_data.fasta',
							'fungi_example_data.fasta',		 
							'plantae_example_data.fasta',
							'protista_example_data.fasta']

	#check outputs were made
	assert sorted(os.listdir('alfie_out')) == fasta_main_outputs
	#tear down the outputs
	for x in fasta_main_outputs:
		os.remove("alfie_out/"+x)
	os.rmdir("alfie_out")

	#process a fastq file
	sys.argv = ['alfie', "-f", ex_fastq_file]
	alf.main()

	fastq_main_outputs = ['animalia_example_data.fastq',
							'bacteria_example_data.fastq',
							'fungi_example_data.fastq',
							'plantae_example_data.fastq',
							'protista_example_data.fastq']

	#check for outputs
	assert sorted(os.listdir('alfie_out')) == fastq_main_outputs

	for x in fastq_main_outputs:
		os.remove("alfie_out/"+x)
	os.rmdir("alfie_out")
