
#import os 

import os
import types
import pytest

from alfie.seqio import file_type, outfile_dict
from alfie.seqio import read_fasta, read_fastq
from alfie.seqio import iter_read_fasta, iter_read_fastq

from alfie import ex_fasta_file, ex_fastq_file


def test_file_type():
	"""Test that the file type is properly identified."""
	assert file_type("file_1.fa") == "fasta"
	assert file_type("file_1.fasta") == "fasta"
	assert file_type("in.file_1.fa") == "fasta"
	assert file_type("file_2.fq") == "fastq"
	assert file_type("file_2.fastq") == "fastq"
	assert file_type("in.file_2.fq") == "fastq"

	with pytest.raises(ValueError):
		file_type("infile_2.txt")

	with pytest.raises(ValueError):
		file_type("in.file_2.csv")


def test_outfile_builder():
	"""Test that the output file set is generated properly."""
	expected_kingdom_dict1 = {0: 'alfie_out/animalia_test.fasta',
						 1: 'alfie_out/bacteria_test.fasta',
						 2: 'alfie_out/fungi_test.fasta',
						 3: 'alfie_out/plantae_test.fasta',
						 4: 'alfie_out/protista_test.fasta'}

	expected_kingdom_dict2 = {0: 'diff_place/animalia_test.fastq',
						 1: 'diff_place/bacteria_test.fastq',
						 2: 'diff_place/fungi_test.fastq',
						 3: 'diff_place/plantae_test.fastq',
						 4: 'diff_place/protista_test.fastq'}

	out1 = outfile_dict("test.fasta")
	assert out1  == expected_kingdom_dict1

	out2 = outfile_dict("in_data/test.fastq", folder_prefix = 'diff_place/') 
	assert out2 == expected_kingdom_dict2

	os.rmdir("alfie_out")
	os.rmdir("diff_place")

def test_fasta_reader():
	""" Test the fasta reader functions."""
	fasta_read = read_fasta(ex_fasta_file)
	
	assert len(fasta_read) == 100
	
	assert fasta_read[0]['name'] == "seq1_plantae"
	assert fasta_read[1]['name'] ==	"seq2_bacteria"
	assert fasta_read[2]['name'] == "seq3_protista"

	assert fasta_read[0]['sequence'][:25] == "TTCTAGGAGCATGTATATCTATGCT"
	assert fasta_read[1]['sequence'][:25] == "ACGGGCTTATCATGGTATTTGGTGC"
	assert fasta_read[2]['sequence'][:25] == "AGTATTAATTCGTATGGAATTAGCA"


def test_fastq_reader():
	""" Test the fastq reader functions."""
	fastq_read = read_fastq(ex_fastq_file)

	assert len(fastq_read) == 100

	for i in range(len(fastq_read)):
		assert list(fastq_read[i].keys()) == ['name', 'sequence', 'strand', 'quality']

	assert fastq_read[0]['sequence'][:25] == "ttctaggagcatgtatatctatgct"
	assert fastq_read[1]['sequence'][:25] == "acgggcttatcatggtatttggtgc"
	assert fastq_read[2]['sequence'][:25] == "agtattaattcgtatggaattagca"


def test_iter_fasta_reader():

	# read in the data
	data = iter_read_fasta(ex_fasta_file, batch = 10)
	
	assert type(data) == types.GeneratorType

	for d in data:
		assert len(d) == 10
		assert list(d[0].keys()) == ['name', 'sequence']


def test_iter_fastq_reader():

	data = iter_read_fastq(ex_fastq_file, batch = 10)

	assert type(data) == types.GeneratorType

	for d in data:
		assert len(d) == 10
		assert list(d[0].keys()) == ['name', 'sequence', 'strand', 'quality']

