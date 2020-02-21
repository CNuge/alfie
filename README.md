# alfie <img src="data/alfie.jpeg" height="140" align="right" alt="Alfie"/>

## alignment free identification of eDNA

## Work in progress

Alfie is a command line tool for kingdom-level classification and processing of DNA sequence data.


## TODO
1. generator read_fasta and read_fastq

2. write_fasta and write_fastq functions
	- take filename and data, append to file

3. decide how to apply the write function, sequentially or in batch?
	- leaning sequentiall



4. Store the 6mer and the 4mer models - have the option to use either from cmd line (make default the 4mer). That way a user can use the more accurate, but slower model if needed