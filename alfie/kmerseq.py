"""
Module containing the KmerFeatures class.

==========
Classes
==========

KmerFeatures : A class to represent a DNA sequence and derive kmer measurements.

"""
import numpy as np

class KmerFeatures:
	"""
	A class to represent a DNA sequence and derive kmer measurements.

	Attributes
	---------
	name : str, the identifier for the sequence.

	k : int, the size of k-mers (substrings of length k) to count. Default is 4.
	
	sequence : str, the nucleotide sequence to generate k-mer counts from. Only counted characters 
		in input are: A, C, G, T. The chatacters N and - are also permitted, but all substrings containing
		these characters are not counted. Presence of any other characters will produce an error.
	
	labels : numpy.ndarray, the different k-mers for the given size of k. Order of k-mers
		is alphabetical and the labels order corresponds to the order of the kmer_freqs values. 
	
	kmer_freqs : numpy.ndarray, the frequencies of the different k-mers. Frequency order
		corresponds to the order of the labels.

	Methods
	---------
	init : takes in a name, sequence, and k value and generates kmer counts and frequency data.

	keys : list, the different k-mers for the given size of k. Order of k-mers
		is alphabetical and the labels order corresponds to the order of the kmer_freqs values. 
	
	values : list, the counts of the different kmers

	freq_values : the frequencies of the different kmers. Calculated as count/total number
		of recorded kmers. Note this means that any substrings with "N" or "-" in them do
		not contribute to the denominator
 	
	items : list, the (key, value) pairs of kmer counts. 
	

	Examples
	---------
	#initiate a class instance
	#by default kmer counts are generated on initialization for k = 4.
	>>> ex_inst = KmerFeatures(name = 'ID1', sequence = 'AAATTTGGGATGGGCCCCACAC')

	# numpy array with the kmer labels
	ex_inst.labels 
	# numpy array with the kmer frequencies
	ex_inst.kmer_freqs

	# count data can be accessed like a dictionary, data returned in list format
	# keys in list format
	>>> ex_inst.keys()
	
	# obtain a list of kmer counts
	>>>	ex_inst.values()

	#or get the values as frequencies
	>>> ex_inst.freq_values()
	
	#obtain the kmer name-count item pairs
	>>> ex_inst.items()

	# You can change the kmer size, and a new set of kmer counts will be 
	# generated
	>>> ex_inst.change_k(2)
	>>> ex_inst.items()
	"""
	def __init__(self, name, sequence, k = 4):

		self.name = name
		self.k = k
		
		up_seq = sequence.upper()
		if self.__check_seq(up_seq) == True:
			self.seq = up_seq

		self.k_dict = self.__kmer_dict(k = self.k)
		self.__count_kmers()

	def __check_seq(self, seq):
		"""Check the input sequence for invalid characters."""
		allowed = {"A", "C", "G", "T", "N", "-"}
		in_set = set(seq)

		if in_set.issubset(allowed) == False:
			raise ValueError("Unallowed characters in input sequence")
		return True		

	def __kmer_build(self, k = 4, dna_list = ["A", "C", "G", "T"]):
		"""Recursive construction of all nucleotide kmer combinations."""
		
		# all the nucleotides to be appended to new kmers
		# note - doing this in alpha order, so keys are pre-sorted
		nts = ['A', 'C', 'G', 'T']
		#decrement k
		k -= 1

		if k == 0:
			return dna_list
		
		else:
			new_dna_list = []

			for kmer in dna_list:
				for n in nts:
					new = kmer + n
					new_dna_list.append(new)
			return self.__kmer_build(k, new_dna_list)

	def __kmer_dict(self, k = 4):
		"""Build an dictionary of all string combos for size k, values all 0. """
		return {k : 0 for k in self.__kmer_build(k)}

	def change_k(self, k, count = True):
		"""Reset k and by default populate the new dictionary. """
		#override existing k
		self.k = k
		#build a new empty dictonary
		self.k_dict = self.__kmer_dict(self.k)

		if count == True:
			self.__count_kmers()

	def __count_kmers(self):
		"""Iterate across a the dna string and count the kmer occurances."""
		for i in range(len(self.seq)-(self.k-1)):
			subseq = self.seq[i:i+self.k]

			if ("N" not in subseq) and ("-" not in subseq):
				self.k_dict[subseq] +=1

	def keys(self):
		"""Returns a list of the kmer keys, in sorted alphabetical order."""
		return [k for k , v in sorted(self.k_dict.items())]

	def values(self):
		"""A list of the kmer count values, maps to keys in alphabetical order."""
		return [v for k , v in sorted(self.k_dict.items())]

	def freq_values(self):
		"""Returns an array of kmer frequencies. """
		vals = np.array(self.values())
		total_count = sum(vals)
		if total_count == 0:
			return vals
		return vals/total_count 

	def items(self):
		return [(k, v) for k, v in sorted(self.k_dict.items())]

	@property
	def labels(self):
		"""The string labels of the kmer frequencies."""
		return np.array(self.keys())

	@property
	def kmer_freqs(self):
		"""A numpy array of kmer frequency values."""
		return np.array(self.freq_values())

