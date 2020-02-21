
import numpy as np

class KmerFeatures:
	"""
	KmerSeq class - obtain kmer counts for a given dna string.
	Only valid characters in input are: A, C, G, T.
	If kmers contain an N or a -, they are skipped over and not counted.
	All other characters will produce an error.

	The default behaviour is to produce an array of kmer frequencies for the
	kmer sizes 6 and 3
	Frequencies are internal to the kmer counts, i.e. the array will sum to 2,
	and the 6mer and 3mer portions of the array both sum to 1

	usage:
	#initiate a class instance
	#by default, kmer counts are generated on initialization.
	ex_inst = KmerFeatures(name = 'ID1',sequence = 'AAATTTGGGATGGGCCCCACAC')

	#obtain the kmer names
	ex_inst.labels #numpy array with the 4160 kmer codings

	#obtain the kmer frequencies
	ex_inst.kmer_freqs

	#You can also set kmer counts for a specific value of k
	ex_inst.get_ks(4) #will build the 4mer dict and populate it
	
	#c
	ex_inst.keys()
	
	#obtain the kmer count values
	ex_inst.values()

	#or the frequencies
	ex_inst.freq_values()
	
	#obtain the kmer name-count item pairs
	ex_inst.items()

	
	"""
	def __init__(self, name, sequence, kmers = [4]):

		self.name = name
		self.kmers = kmers
		self.seq = sequence.upper()

		self.outkeys = []
		self.outvals = []

		for k in kmers:
			self.k = k
			self.k_dict = self.kmer_dict(self.k)
			self.count_kmers()

			self.outkeys.extend(self.keys())
			self.outvals.extend(self.freq_values())


	def kmer_build(self, k = 6, dna_list = ['A', 'C', 'G', 'T']):
		"""
		recursive construction of all nucleotide kmer combinations 
		"""
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
			return self.kmer_build(k, new_dna_list)


	def kmer_dict(self, k = 6):
		"""
		builds a dictonary where the keys are all the nt combinations
		for the specified k and the values are integers
		"""
		return {k : 0 for k in self.kmer_build(k)}


	def get_ks(self, k, count = True):
		"""
		user sets or resets k, and the corresponding dict for 
		the given class instance. by default the dict is populated
		"""
		#override existing k
		self.k = k
		#build a new empty dictonary
		self.k_dict = self.kmer_dict(self.k)

		if count == True:
			self.count_kmers()


	def count_kmers(self):
		"""
		iterate across a the dna string and count the kmer occurances
		this function returns nothing, just populates the kmer dictonary
		"""
		for i in range(len(self.seq)-(self.k-1)):
			subseq = self.seq[i:i+self.k]

			if ("N" not in subseq) and ("-" not in subseq):
				self.k_dict[subseq] +=1

	def keys(self):
		""" 
		returns a list of the kmer keys, in sorted alphabetical order
		"""

		return [k for k , v in sorted(self.k_dict.items())]


	def values(self):
		"""
		returns a list of the kmer values, corresponding to the keys in sorted 
		alphabetical order
		"""
		return [v for k , v in sorted(self.k_dict.items())]

	def freq_values(self):
		vals = np.array(self.values())
		total_count = sum(vals)
		if total_count == 0:
			return vals
		return vals/total_count 

	def items(self):
		return [(k, v) for k, v in sorted(self.k_dict.items())]

	@property
	def labels(self):
		""" the labels of the kmer freq columns """
		return np.array(self.outkeys)

	@property
	def kmer_freqs(self):
		""" the kmer frequency values """

		return np.array(self.outvals)
