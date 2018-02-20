'''
Automate analysis of Sanger trace data.  


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

'''


from sys import argv
from Bio import SeqIO
import os
from subprocess import call


# Parameters

default_blast_db = '/ebio/abt3_projects/databases/NCBI_blastdb/nt'
min_qual_cutoff    = 50
max_qual_cutoff      = 150
min_region_length_cutoff = 550


# args/set-up

input_abi_file = os.path.abspath(argv[1])
if input_abi_file[-3:] != 'ab1':
	raise ValueError('Input must be .ab1 file')
output_prefix = input_abi_file[:-4]
sample_name = output_prefix.split('/')[-1]

call('mkdir '+sample_name,shell=True)
output_prefix = '%s/%s/%s' % (os.getcwd(),sample_name,sample_name)
print(output_prefix)

mode = argv[2]
if mode not in ('known','unknown'):
	raise ValueError('Must select mode: either "known" or "unknown"')

if mode == 'known':
	sequence_to_align = argv[3]

if mode == 'unknown':
	blast_db = argv[3]
	if argv[3].lower() == 'default':
		blast_db = default_blast_db

verbose = argv[4]
if verbose in ('V','v'):
	verbose = True
elif verbose not in ('S','s'):
	raise ValueError("Must select either 'V'erbose or 'S'hort for output")



# Opens the abi file and grabs the phred scores, and the sequence itself

abi_file = list(SeqIO.parse(open(input_abi_file,"rb"), "abi"))[0]
phred_quality = abi_file.letter_annotations["phred_quality"]
seq = abi_file.seq



# Grab longest continuous stretch of good quality bp


for qual_range in range(max_qual_cutoff,min_qual_cutoff,-1):
	good_region = False
	region_start = 0
	good_qual_regions = []
	for i in range(len(seq)-3):
		if phred_quality[i]+phred_quality[i+1]+phred_quality[i+2] < qual_range:
			if good_region == True:
				region_end = i
				good_region_DNA = seq[region_start:region_end]
				good_region_length = len(good_region_DNA)
				good_region_phred_sum = sum(phred_quality[region_start:region_end])
				avg_phred = good_region_phred_sum*1.0 / good_region_length
				good_qual_regions.append([good_region_length, str(good_region_DNA), avg_phred])
				good_region = False
		if phred_quality[i]+phred_quality[i+1]+phred_quality[i+2] >= qual_range:
			if good_region == False:
				region_start = i
				good_region = True
	sorted_good_qual_regions = sorted(good_qual_regions, reverse=True)
	if not sorted_good_qual_regions[0]:
		print('No good quality region')
		sys.exit()
	longest_good_qual = sorted_good_qual_regions[0]
	if longest_good_qual[0] >= min_region_length_cutoff:
		print('\nLongest high quality stretch was',longest_good_qual[0],' bp, and had an avg Phred of',str(round(longest_good_qual[2],2)),'\n')
		break

tmp_sanger = output_prefix+'_tmp_sanger.fasta'
tmp_sanger_output = open(tmp_sanger,"w")
tmp_sanger_output.write('>%s\n%s\n' % (sample_name,longest_good_qual[1]))
tmp_sanger_output.close()





# Known mode

if mode == 'known':

	tmp_to_align = output_prefix+'_tmp_to_align.fasta'
	tmp_to_align_output = open(tmp_to_align,'w')
	tmp_to_align_output.write('>Parent_sequence\n'+str(sequence_to_align))
	tmp_to_align_output.close()
	

# Performs MegaBlast with Sanger read and input sequence

	call('blastn -query '+tmp_sanger+' -subject '+tmp_to_align+' -outfmt 6 -out '+output_prefix+'_BlastResults.txt 2>/dev/null',shell=True)
	for l in open(output_prefix+'_BlastResults.txt','U'):
		E_value = l.split("\t")[10]

# Performs MAFFT alignment of input and Sanger sequence

	tmp_mafft = output_prefix+'_tmp_mafft.fasta'
	tmp_mafft_output = open(tmp_mafft,'w')
	tmp_mafft_output.write('\n'.join(('>'+sample_name,str(longest_good_qual[1]),'>Parent',str(sequence_to_align))))
	tmp_mafft_output.close()

	mafft_file_name = output_prefix+'_alignment.txt'
	call('mafft --clustalout '+tmp_mafft+' > '+mafft_file_name+' 2>/dev/null',shell=True)

	
	if verbose == True:
		for l in open(mafft_file_name,'U'):
			print(l)
	print('\n %s blastn E-value to input seq was %s\n' % (sample_name,str(E_value)))


# Clean-up temp files and move results to correct folder

	call('rm %s %s' % (tmp_mafft,tmp_to_align,shell=True)


# Unknown Mode

if mode == 'unknown':
	call('blastn -query '+tmp_sanger+' -db '+blast_db+' -num_alignments 50 -num_descriptions 50 -outfmt 0 -out '+output_prefix+'_BlastResults.txt 2>/dev/null',shell=True)
	
	if verbose == True:
		print_now = False
		for l in open(output_prefix+'_BlastResults.txt','U'):
			if l[0] == '>' and print_now == True:
				break
			if l[0] == '>':
				print_now = True
			if print_now == True:
				print(l)
