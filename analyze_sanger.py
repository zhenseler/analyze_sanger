'''
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

'''


import sys
from Bio import SeqIO
import os
from subprocess import call


# Parameters

default_blast_db = '/ebio/abt3_projects/databases/NCBI_blastdb/nt'
min_qual_cutoff    = 10
max_qual_cutoff      = 150
min_region_length_cutoff = 550


# args/set-up

input_abi_file = os.path.abspath(sys.argv[1])
if input_abi_file[-3:] != 'ab1':
	raise ValueError('Input must be .ab1 file')
output_prefix = input_abi_file[:-4]
sample_name = output_prefix.split('/')[-1]

if os.path.isdir(sample_name) == False:
	call('mkdir ' + sample_name, shell = True)
output_prefix = '%s/%s/%s' % (os.getcwd(), sample_name, sample_name)

mode = sys.argv[2]
if mode not in ('known', 'unknown_local', 'unknown_remote'):
	raise ValueError('Must select mode: either "known" "unknown_local" or "unknown_remote"')

if mode == 'known':
	sequence_to_align = sys.argv[3]

if mode == 'unknown_local':
	blast_db = sys.argv[3]
	if sys.argv[3].lower() == 'default':
		blast_db = default_blast_db

if mode == 'unknown_remote':
	blast_db = sys.argv[3]
	if sys.argv[3].lower() == 'default':
		blast_db = 'nt' # By default use nt when doing unknown


# Opens the abi file and grabs the phred scores, and the sequence itself

abi_file = list(SeqIO.parse(open(input_abi_file,"rb"), "abi"))[0]
phred_quality = abi_file.letter_annotations["phred_quality"]
seq = abi_file.seq

if seq == "NNNNN":
	print(sample_name + "\tNNNNN\tNA\tNA\tNA\n")
	sys.exit()



# Grab longest continuous stretch of good quality bp

# Start at highest quality cutoff, and decrease to minimum if needed
for qual_range in range(max_qual_cutoff, min_qual_cutoff,-1):
	good_region = False
	region_start = 0
	good_qual_regions = []
	
	# Scan through sequence
	for i in range(len(seq)-3):
		# If in a good stretch, and quality dips below cutoff
		if phred_quality[i] + phred_quality[i+1] + phred_quality[i+2] < qual_range:
			# If there was any stretch that wasn't trash
			if good_region == True:
				region_end = i
				good_region_DNA = seq[region_start:region_end]
				good_region_length = len(good_region_DNA)
				good_region_phred_sum = sum(phred_quality[region_start:region_end])
				avg_phred = good_region_phred_sum * 1.0 / good_region_length
				good_qual_regions.append([good_region_length, str(good_region_DNA), avg_phred])
				good_region = False
		# If quality is above cutoff
		elif phred_quality[i] + phred_quality[i+1] + phred_quality[i+2] >= qual_range:
			if good_region == False:
				region_start = i
				good_region = True
				
sorted_good_qual_regions = sorted(good_qual_regions, reverse = True)
if sorted_good_qual_regions == [] or sorted_good_qual_regions[0][0] < min_region_length_cutoff:
	print(sample_name + "\tBelow quality threshold\tNA\tNA\tNA\n")
	sys.exit()
longest_good_qual = sorted_good_qual_regions[0]
if longest_good_qual[0] >= min_region_length_cutoff:
	final_average_phred = str(round(longest_good_qual[2], 2))

tmp_sanger = output_prefix + '_tmp_sanger.fasta'
tmp_sanger_output = open(tmp_sanger,"w")
tmp_sanger_output.write('>%s\n%s\n' % (sample_name, longest_good_qual[1]))
tmp_sanger_output.close()





# Known mode

if mode == 'known':

	blast_result = ""

	tmp_to_align = output_prefix +' _tmp_to_align.fasta'
	tmp_to_align_output = open(tmp_to_align,'w')
	tmp_to_align_output.write('>Parent_sequence\n' + str(sequence_to_align))
	tmp_to_align_output.close()
	

# Performs MegaBlast with Sanger read and input sequence

	blast_output_name = output_prefix + '_blast_results.txt'
	call('blastn -query ' + tmp_sanger +' -subject ' + tmp_to_align+' -outfmt 6 -out '\
			+ blast_output_name + ' 2>/dev/null', shell = True)
	if os.path.getsize(blast_output_name) <= 0:
		blast_result = "Below default cutoff"
	else:
		E_value_list = []
		for l in open(blast_output_name,'U'):
			E_value_list.append(l.split("\t")[10])
		blast_result = " ,".join(E_value_list)


# Performs MAFFT alignment of input and Sanger sequence

	tmp_mafft = output_prefix + '_tmp_mafft.fasta'
	tmp_mafft_output = open(tmp_mafft,'w')
	tmp_mafft_output.write('\n'.join(('>'+sample_name, longest_good_qual[1], \
							'>Parent', sequence_to_align)))
	tmp_mafft_output.close()

	mafft_file_name = output_prefix + '_alignment.txt'
	call('mafft --clustalout ' + tmp_mafft+' > ' + mafft_file_name + \
			' 2>/dev/null',shell=True)



# Clean-up temp files and move results to correct folder

	call('rm %s %s' % (tmp_mafft,tmp_to_align), shell = True)

	print("\t".join([sample_name, longest_good_qual[1], str(longest_good_qual[0]), \
			final_average_phred, blast_result]) + "\n")


# Unknown modes

if mode == 'unknown_local':

	blast_output_name = output_prefix + '_blast_results.txt'
	
	call('blastn -query ' + tmp_sanger + ' -db ' + blast_db + \
			' -num_alignments 50 -num_descriptions 50 -outfmt 0 -out '\
			 + blast_output_name + ' 2>/dev/null', shell = True)
			 
	for l in open(blast_output_name, 'U'):
		if l[0] == '>':
			print("\t".join([sample_name, longest_good_qual[1], str(longest_good_qual[0]), \
				final_average_phred, l]) + '\n')
			
				
if mode == 'unknown_remote':
	blast_output_name = output_prefix + '_blast_results.txt'
	
	call('blastn -query ' + tmp_sanger + ' -db ' +  blast_db + ' -remote \
			-outfmt 6 -max_target_seqs 5 -out ' + blast_output_name, shell = True)
	
	description_list = []
	tmp_esearch = "tmp_esearch.fa"
	for l in open(blast_output_name, 'U'):
		locus = l.split("\t")[1]
		call('esearch -db nucleotide -query ' + locus + ' | efetch -format fasta > ' + tmp_esearch, shell = True) 
		for esearch_l in open(tmp_esearch, 'U'):
			if(esearch_l.startswith('>')):
				description_list.append(esearch_l.strip())
	
	print("\t".join([sample_name,  '\t'.join(description_list), str(longest_good_qual[0]), final_average_phred, longest_good_qual[1]]))