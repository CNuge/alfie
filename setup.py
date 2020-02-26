from setuptools import setup, find_packages

with open('requirements.txt') as f:
	requirements = f.readlines()
	requirements = [x.rstrip() for x in requirements]

long_description = "Alfie is a package and CLI program, what is does will be written here."

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


^ all the above works, but it is not storing the model in the executable, 
need a way to make this get installed into the package.
"""