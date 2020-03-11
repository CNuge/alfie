from setuptools import setup, find_packages

with open('requirements.txt') as f:
	requirements = f.readlines()
	requirements = [x.rstrip() for x in requirements]

long_description = """
alfie: an alignment-free, kingdom level taxonomic classifier for DNA barcode data. 

Alfie classifies sequences using a neural network which takes k-mer frequencies (default k = 4)
as inputs and makes kingdom level classification predictions. At present, the program contains 
trained models for classification of cytochrome c oxidase I (COI) barcode sequences to the 
taxonomic level: kingdom. The program is effective at classifying sequences >200 base pairs in 
length, and no alignment information is needed.

Alfie can be deployed from the command line for rapid file-to-file classification of sequences. 
This is an effective means of separating contaminant sequences in a DNA metabarcoding or 
environmental DNA dataset from sequences of interest.

For increased control, alfie can also be deployed as a module from within Python. The alfie 
module contains functions that can aid a user in the training and application of a custom 
alignment-free classifier, which allows the program to be applied to different DNA barcodes 
(or genes) or on different taxonomic levels.
"""


setup(
	name = 'alfie',
	version = '0.1',
	author = 'Cam Nugent',
	author_email = 'nugentc@uoguelph.ca',
	url = 'https://github.com/CNuge/alfie',
	description = 'alignment free identification of edna',
	long_description = long_description,
	license= 'LICENSE.md',
	packages = find_packages(),
	package_data={'alfie': ['data/*']},
	entry_points = {
	'console_scripts':[
	'alfie = alfie.alf:main']
	},
	python_requires='>=3.6',
	install_requires = requirements,

	)


"""
create the release:
python setup.py sdist
install the release:
python3 setup.py install

#then check from home dir if the package works with
alfie -h

#can check to see if functions available with:
from alfie.kmerseq import KmerFeatures
?KmerFeatures
"""