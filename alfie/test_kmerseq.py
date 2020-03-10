
import unittest

from kmerseq import KmerFeatures

class KmerTests(unittest.TestCase):
	"""Unit tests for the KmerFeatures class."""
	@classmethod
	def setUpClass(self):
		"""Initiate the test class instance."""
		self.test_kmers = KmerFeatures("test1", 
							"aaaaaattttttatatatgcgcgccccccgccgcgccgggc")

	def test_KmerFeatures(self):
		
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

		with self.assertRaises(ValueError):
			self.assertEqual(KmerFeatures("test1", "NOTDNA"))


if __name__ == '__main__':
	unittest.main()

