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
warnings.filterwarnings("ignore", category=RuntimeWarning) 

from Functions import CustomRemovalSoftClippedBases2

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
		print ('There is a Results folder for the mature genome alignment already created for this sample ! Maybe the analysis has already been done (Cheek the Results folder)')
else:
	print ('The results folder for this sample has not been created so the first step has not been done ! ')
	raise SystemExit

os.chdir('../Results/'+sample_name+'/Alignment_PG/')

# Alignment
print ('3- Performing the alignment against the Precursor tRNA Genome')
command = "bowtie2 --local -p " + str(nucleos)
os.system(command + ' -N 0 -x '+refgen_path+'precursor_fam_tRNA_refgenome ../Alignment_MG/'+sample_name+'_WGloc_only_trna_precursor_and_MGunmapped.fastq'+' --un '+sample_name+'_unmapped_PGloc.fastq'+' | samtools view -bSF4 - > '+sample_name+'_PGloc_mapped.bam')

# Processing files (sorting and indexing).
os.system('samtools sort '+sample_name+'_PGloc_mapped.bam'+ ' -o ' +sample_name+'_PGloc_mapped_sort.bam')
os.system('samtools index '+sample_name+'_PGloc_mapped_sort.bam') 

# Removing soft clipped bases
CustomRemovalSoftClippedBases2(sample_name,'_PGloc_mapped_')

# Processing files (sorting and indexing).
os.system('samtools sort '+sample_name+'_PGloc_mapped_soft_clipped_rm.bam'+ ' -o ' +sample_name+'_PGloc_mapped_soft_clipped_rm_sort.bam')
os.system('samtools index '+sample_name+'_PGloc_mapped_soft_clipped_rm_sort.bam') 


#Add the file to the final result folder
if not os.path.exists('../Final_results'):
	os.mkdir('../Final_results')
os.system('cp -R '+sample_name+'_PGloc_mapped_soft_clipped_rm.bam ../Final_results/') 
os.rename('../Final_results/'+sample_name+'_PGloc_mapped_soft_clipped_rm.bam', '../Final_results/'+sample_name+'_PGloc_mapped.bam')

# Processing files (sorting and indexing).
os.system('samtools sort '+'../Final_results/'+sample_name+'_PGloc_mapped.bam'+ ' -o ' +'../Final_results/'+sample_name+'_PGloc_mapped_sort.bam')
os.system('samtools index '+'../Final_results/'+sample_name+'_PGloc_mapped_sort.bam')
sys.exit()
