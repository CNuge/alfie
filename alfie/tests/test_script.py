import pytest

from alfie import alf

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


