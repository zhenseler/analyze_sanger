# analyze_sanger
Automate analysis of Sanger trace data


Dependencies:

	- blast
	- Biopython
	- mafft

Usage:

	- python analyze_sanger.py <input_file> <mode> <parent/db> <verbosity>
	
		- input_file: 	input .ab1 Sanger trace file
		
		- mode: 		'known': 	compare Sanger to input sequence
						'unknown': 	blast Sanger against database
						
		- parent/db: 	if mode = 'known', DNA sequence string to compare Sanger to.
						if mode = 'unknown', path to blast database
						
		- verbosity:	'v': verbose, prints out alignment
						's': short, no print statements
	
Workflow (known):
	
	1.	Grab largest stretch of trace above quality cutoff
	2.	blastn against input parent sequence. Print e-value
	3.	Align to input parent sequence via MAFFT
			Print alignment if verbosity = 'v'
	
Workflow (unknown):

	1. Grab largest stretch of trace above quality cutoff
	2. blastn against input blast database
			Print alignment of top if if verbosity = 'v'

To-do:

	- 'no output' mode: no output files created
	- 'fasta convert' mode: just does fasta conversion
