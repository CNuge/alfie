import unittest

from seqio import file_type


class SeqioTests(unittest.TestCase):
	"""
	unit tests for the io functions associated with the main alfie executable
	"""
	def test_file_type(self):
		self.assertEqual(file_type("file_1.fa") , "fasta")
		self.assertEqual(file_type("file_1.fasta") , "fasta")
		self.assertEqual(file_type("in.file_1.fa") , "fasta")
		self.assertEqual(file_type("file_2.fq") , "fastq")
		self.assertEqual(file_type("file_2.fastq") , "fastq")
		self.assertEqual(file_type("in.file_2.fq") , "fastq")

		with self.assertRaises(ValueError):
			self.assertEqual(file_type("infile_2.txt"))
	
		with self.assertRaises(ValueError):
			self.assertEqual(file_type("in.file_2.csv"))





if __name__ == '__main__':
	unittest.main()