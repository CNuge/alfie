import pytest
from alfie.kmerseq import KmerFeatures

def test_KmerFeatures():
	"""Unit tests for the KmerFeatures class."""
	test_kmers = KmerFeatures("test1", 
						"aaaaaattttttatatatgcgcgccccccgccgcgccgggc")
	
	assert test_kmers.name == "test1"

	assert test_kmers.labels.shape == (256,)

	assert list(test_kmers.labels[:3]) == ["AAAA", "AAAC", "AAAG"]

	assert list(test_kmers.labels[-3:]) == ["TTTC", "TTTG", "TTTT"]

	assert test_kmers.kmer_freqs.shape == (256,)

	assert test_kmers.items()[0] == ("AAAA", 3)

	with pytest.raises(ValueError):
		KmerFeatures("test1", "NOTDNA")


	#chage the k, asser it was done properly
	test_kmers.change_k(2)

	assert test_kmers.kmer_freqs.shape == (16,)
	assert list(test_kmers.labels[:4]) == ["AA", "AC", "AG", "AT"]
	assert list(test_kmers.labels[-4:]) == ["TA", "TC", "TG", "TT"]


	#test handling of empty kmer counts
	#should make all the frequencies 0, while avoiding a divide by zero.
	e_test = KmerFeatures("test2", "")

	assert list(e_test.freq_values()) == [0]*256