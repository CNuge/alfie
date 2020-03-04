
import os 
import unittest

from seqio import file_type, outfile_dict, read_fasta, read_fastq


class SeqioTests(unittest.TestCase):
	"""
	unit tests for the io functions associated with the main alfie executable
	"""
	@classmethod
	def setUpClass(self):
		"""
		initiate the test class instance with the 
		"""
		self._expected_kingdom_dict = {0: 'alfie_out/animalia_test.fasta',
										 1: 'alfie_out/bacteria_test.fasta',
										 2: 'alfie_out/fungi_test.fasta',
										 3: 'alfie_out/plantae_test.fasta',
										 4: 'alfie_out/protista_test.fasta'}
		
		self._fasta_infile = 'data/small_unittest.fasta'
		self._fastq_infile = 'data/small_unittest.fastq'

	@classmethod
	def tearDown(self):
		"""
		after unit tests, remove the temporary outputs
		"""
		try:
			os.rmdir("alfie_out")
		except OSError:
			pass

	def test_file_type(self):
		"""
		test that the file type is properly identified
		"""
		self.assertEqual(file_type("file_1.fa"), 
						"fasta")
		self.assertEqual(file_type("file_1.fasta"),
						"fasta")
		self.assertEqual(file_type("in.file_1.fa"),
						"fasta")
		self.assertEqual(file_type("file_2.fq"),
						"fastq")
		self.assertEqual(file_type("file_2.fastq"),
						"fastq")
		self.assertEqual(file_type("in.file_2.fq"),
						"fastq")

		with self.assertRaises(ValueError):
			self.assertEqual(file_type("infile_2.txt"))
	
		with self.assertRaises(ValueError):
			self.assertEqual(file_type("in.file_2.csv"))


	def test_outfile_builder(self):
		"""
		test that the output file set is generated properly
		"""
		self.assertEqual(outfile_dict("test.fasta"), 
						self._expected_kingdom_dict)

		self.assertEqual(outfile_dict("in_data/test.fasta"), 
				self._expected_kingdom_dict)


	def test_fasta_reader(self):
		self._fasta_read = read_fasta(self._fasta_infile)
		
		self.assertEqual(len(self._fasta_read), 3)
		
		self.assertEqual(self._fasta_read[0]['name'], 
						"example_read_1")
		self.assertEqual(self._fasta_read[1]['name'], 
						"example_read_2")
		self.assertEqual(self._fasta_read[2]['name'], 
						"example_read_3")

		self.assertEqual(self._fasta_read[0]['sequence'][:25], 
						"gttcaacaaatcataaagatattgg")
		self.assertEqual(self._fasta_read[1]['sequence'][:25],
						"attcaaccaatcataaagatattgg")
		self.assertEqual(self._fasta_read[2]['sequence'][:25],
						"gatatagcatttccacgacttaata")

	def test_fastq_reader(self):
		self._fastq_read = read_fastq(self._fastq_infile)

		self.assertEqual(len(self._fastq_read), 3)
	
		for i in range(len(self._fastq_read)):
			self.assertEqual(list(self._fastq_read[i].keys()),
							['name', 'sequence', 'strand', 'quality'])

		self.assertEqual(self._fastq_read[0]['sequence'][:25], 
						"gttcaacaaatcataaagatattgg")
		self.assertEqual(self._fastq_read[1]['sequence'][:25],
						"attcaaccaatcataaagatattgg")
		self.assertEqual(self._fastq_read[2]['sequence'][:25],
						"cggaattcccggatcaataattggg")


if __name__ == '__main__':
	unittest.main()

