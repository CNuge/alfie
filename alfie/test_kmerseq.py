
import unittest

from kmerseq import KmerFeatures

class KmerTests(unittest.TestCase):
	"""
	unit tests for the io functions associated with the main alfie executable
	"""
	@classmethod
	def setUpClass(self):
		"""
		initiate the test class instance with the 
		"""
		self.test_kmers = KmerFeatures("test1", 
							"aaaaaattttttatatatgcgcgccccccgccgcgccgggc")

	def test_file_type(self):
		"""
		test that the file type is properly identified
		"""
		self.assertEqual(self.test_kmers.name, 
						"test1")

		self.assertEqual(self.test_kmers.labels.shape,
						(256,))

		self.assertEqual(list(self.test_kmers.labels[:3]),
				['AAAA', 'AAAC', 'AAAG'])

		self.assertEqual(list(self.test_kmers.labels[-3:]),
				['TTTC', 'TTTG', 'TTTT'])

		self.assertEqual(self.test_kmers.kmer_freqs.shape,
						(256,))


if __name__ == '__main__':
	unittest.main()

