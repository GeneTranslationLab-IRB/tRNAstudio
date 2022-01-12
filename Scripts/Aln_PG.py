#!/usr/bin/env python3

'''
This script performs the alignment against the precursor genome.
'''

###############################################################################
# Built-in/Generic Imports 
import os 
import sys
import subprocess
import multiprocessing
import warnings
import Bio
from Bio import SeqIO
import pysam
warnings.filterwarnings("ignore", category=RuntimeWarning) 

from Functions import custom_removal_soft_clipped_bases_2

nucleos = multiprocessing.cpu_count()
###############################################################################

# Arguments
sample_name=(sys.argv[1].split('.'))[0]
path=(os.path.dirname(os.path.realpath(__file__))).split('/')[:-1]
refgen_path=('/').join(path)+'/Reference_Genomes/'

if os.path.exists('../Results/'+sample_name):
	if not os.path.exists('../Results/'+sample_name+'/Alignment_PG'):
			os.makedirs('../Results/'+sample_name+'/Alignment_PG')
	else:
		print ('There is a Results folder for the mature genome alignment already created for this sample ! Maybe the analysis has already been done (Cheek the Results folder)! Reanalyzing the data')
else:
	print ('The results folder for this sample has not been created so the first step has not been done! ')
	raise SystemExit

os.chdir('../Results/'+sample_name+'/Alignment_PG/')

# Alignment
print ('3- Performing the alignment against the Precursor tRNA Genome')
command = "bowtie2 --local -p " + str(nucleos)
os.system(command + ' -N 0 -x '+refgen_path+'precursor_fam_tRNA_refgenome ../Alignment_MG/'+sample_name+'_WGloc_only_trna_precursor_and_MGunmapped.fastq'+' --un '+sample_name+'_unmapped_PGloc.fastq'+' | samtools view -bSF4 - > '+sample_name+'_precursor.bam')

# Processing files (sorting and indexing)
os.system('samtools sort '+sample_name+'_precursor.bam'+ ' -o ' +sample_name+'_precursor_sort.bam')
os.system('samtools index '+sample_name+'_precursor_sort.bam') 

# Filtering reads that are not mapping precursor regions
precursor_seq = SeqIO.parse(refgen_path+'precursor_fam_tRNA_refgenome.fa', "fasta")
bamfile = pysam.AlignmentFile(sample_name+'_precursor_sort.bam', "rb")
out_precursor_filtered = pysam.AlignmentFile(sample_name+'_PGloc_mapped_filtered.bam', "wb", template=bamfile)
out_mature_filtered = pysam.AlignmentFile(sample_name+'_mature_mapped_filtered.bam', "wb", template=bamfile)

for tRNA in precursor_seq:
	tRNA_len = len(tRNA.seq)
	for read in bamfile.fetch(tRNA.id):
		readStart, readEnd = read.reference_start, read.reference_end
		if readStart > 49 and readEnd < tRNA_len-49:
			out_mature_filtered.write(read)
		else:
			out_precursor_filtered.write(read)
out_precursor_filtered.close()
out_mature_filtered.close()
bamfile.close()

os.system('samtools sort '+sample_name+'_PGloc_mapped_filtered.bam'+ ' -o ' +sample_name+'_PGloc_mapped_filtered_sort.bam')
os.system('samtools index '+sample_name+'_PGloc_mapped_filtered_sort.bam') 

# Removing soft clipped bases
custom_removal_soft_clipped_bases_2(sample_name,'_PGloc_mapped_filtered_')

# Processing files (sorting and indexing)
os.system('samtools sort '+sample_name+'_PGloc_mapped_filtered_soft_clipped_rm.bam'+ ' -o ' +sample_name+'_PGloc_mapped_filtered_soft_clipped_rm_sort.bam')
os.system('samtools index '+sample_name+'_PGloc_mapped_filtered_soft_clipped_rm_sort.bam') 

# Adding the file to the final result folder
os.system('cp -R '+sample_name+'_PGloc_mapped_filtered_soft_clipped_rm_sort.bam ../Alignments/') 
os.rename('../Alignments/'+sample_name+'_PGloc_mapped_filtered_soft_clipped_rm_sort.bam', '../Alignments/'+sample_name+'_precursor.bam')

# Processing files (sorting and indexing)
os.system('samtools sort '+'../Alignments/'+sample_name+'_precursor.bam'+ ' -o ' +'../Alignments/'+sample_name+'_precursor_sort.bam')
os.system('samtools index '+'../Alignments/'+sample_name+'_precursor_sort.bam')


os.system('samtools fastq '+sample_name+'_mature_mapped_filtered.bam > '+sample_name+'_mature_mapped_filtered.fastq')
os.system('cat '+sample_name+'_mature_mapped_filtered.fastq '+sample_name+'_unmapped_PGloc.fastq > '+sample_name+'_mature_and_PG_unmapped.fastq')

print ('\n')
