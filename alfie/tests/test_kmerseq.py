import pytest
from alfie.kmerseq import KmerFeatures

def test_KmerFeatures():
	"""Unit tests for the KmerFeatures class."""
	test_kmers = KmerFeatures("test1", 
						"aaaaaattttttatatatgcgcgccccccgccgcgccgggc")
	
	assert test_kmers.name == "test1"

	assert test_kmers.labels.shape == (256,)

	assert list(test_kmers.labels[:3]) == ['AAAA', 'AAAC', 'AAAG']

	assert list(test_kmers.labels[-3:]) == ['TTTC', 'TTTG', 'TTTT']

	assert test_kmers.kmer_freqs.shape == (256,)

	with pytest.raises(ValueError):
		KmerFeatures("test1", "NOTDNA")



