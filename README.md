# analyze_sanger
Automate analysis of Sanger trace data.  


Dependencies:

	- blast				(conda install -c bioconda blast)
	- Biopython 		(conda install -c conda-forge biopython)
	- mafft 			(conda install -c bioconda mafft)
	- entrez-direct		(conda install -c bioconda entrez-direct)		

Usage:

	- python analyze_sanger.py <input_file> <mode> <parent/db>
	
		- input_file: 	input .ab1 Sanger trace file
		
		- mode: 		'known': 	compare Sanger to input sequence
						'unknown_local': 	blast Sanger against local database
						'unknown_remote': blast Sanger against remote database
						
		- parent/db: 	if mode = 'known', DNA sequence string to compare Sanger to.
						if mode = 'unknown_local', path to local blast database
						if mode = 'unknown_remote', which remote db to use

						
		** Want to analyze multiple files and put the output into a single file? **
		
		   Use bash and pipe printed output to file
		   for file in <folder>; do python <script location> $file <mode> 
		   		<sequence to match to>; done > <desired output file>
	
	
Workflow (known):
	
	1.	Descending from max quality cutoff to input minimum, grab longest stretch of 
			DNA above that quality (see Parameters below). Stop when a sequence is 
			obtained longer than input minimum sequence length (see Parameters)
	2.	blastn against input parent sequence with default parameters
	3.	Align to input parent sequence via MAFFT
			
	
Workflow (unknown):

	1. Grab largest stretch of trace above quality cutoff
	2. blastn against input blast database (either local or remote)

To-do:

	- 'no output' mode: no output files created
	- 'fasta convert' mode: just does fasta conversion

To-do:

	- 'no output' mode: no output files created
	- 'fasta convert' mode: just does fasta conversion
