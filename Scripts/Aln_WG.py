#!/usr/bin/env python3

'''
This script performs the first step of the pipeline 
- Alignment against the whole genome
- Select reads mapping to tRNA genes
- Removal of soft clipped bases 
- Classification between precursor and processed/mature reads.
'''

###############################################################################
# Built-in/Generic Imports and Modules #
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning) 

import argparse
import subprocess
import os 
import sys
import glob
import re
import multiprocessing
import pysam
import pandas as pd
from Bio import SeqIO

# Import functions #
from Functions import *

nucleos = multiprocessing.cpu_count()

###############################################################################

#####  All the required arguments  ######

# Samples path
samples_path = sys.argv[1]
samples = os.listdir(samples_path)
sample = sys.argv[2]
sample_name = sample.split('.')[0]

print ('ANALYZING SAMPLE: '+sample)

# Ref genomes path
refgen_path = '../Reference_Genomes/'

# tRNA Coordinates
tRNA_coordinates = pd.read_csv(refgen_path+'info/tRNAsCoordinates.bed', sep='\t', header=None)
header = ['chrom', 'Start', 'End', 'name', 'score', 'strand']
tRNA_coordinates.columns = header[:len(tRNA_coordinates.columns)]

# Intron Coordinates
intron_coordinates = pd.read_csv(refgen_path+'info/intron.bed', sep='\t', header=None)
header = ['chrom', 'Start', 'End', 'name', 'strand']
intron_coordinates.columns = header[:len(intron_coordinates.columns)]

# Create sample folder in Results folder
if not os.path.exists('../Results/'+sample_name+'/Alignment_WG/'):
            os.makedirs('../Results/'+sample_name+'/Alignment_WG/')  

if not os.path.exists('../Results/'+sample_name+'/Alignments/'):
            os.makedirs('../Results/'+sample_name+'/Alignments/')    

# Bowtie2 command
command = "bowtie2 --local -p "+ str(nucleos)

# Alignment against the Whole Genome with Bowtie2
print ('1- Performing the alignment against the Whole Human Genome:')
os.system(command + ' -N 0 -x '+refgen_path+'genome '+samples_path+ '/'+sample+' --un ../Results/'+sample_name+'/Alignment_WG/'+sample_name+'_unmapped_WGloc.fastq | samtools view -bSF4 - > ../Results/'+sample_name+'/Alignment_WG/'+sample_name+'_WGloc_mapped.bam')

# Pipeline to classify between precursor and processed reads. 
# Set working directory.
os.chdir('../Results/'+sample_name+'/Alignment_WG/')

# Processing files (sorting and indexing)
os.system('samtools sort '+sample_name+'_WGloc_mapped.bam'+ ' -o ' +sample_name+'_WGloc_mapped_sort.bam')
os.system('samtools index '+sample_name+'_WGloc_mapped_sort.bam')   

print ('Selecting reads mapping to tRNAs...')
os.system('samtools view -b -L '+'../../../Reference_Genomes/info/tRNAsCoordinates.bed'+' '+sample_name+'_WGloc_mapped.bam > '+sample_name+'_WGloc_only_trna.bam')
os.system('samtools sort '+sample_name+'_WGloc_only_trna.bam'+ ' -o ' +sample_name+'_WGloc_only_trna_sort.bam')
os.system('samtools index '+sample_name+'_WGloc_only_trna_sort.bam') 

print ('Custom removal of soft clipped bases...')
custom_removal_soft_clipped_bases(sample_name,'_WGloc_only_trna_', tRNA_coordinates) 
os.system('samtools sort '+sample_name+'_WGloc_only_trna_soft_clipped_rm_all_chr.bam'+ ' -o ' +sample_name+'_WGloc_only_trna_soft_clipped_rm_all_chr_sort.bam')
os.system('samtools index '+sample_name+'_WGloc_only_trna_soft_clipped_rm_all_chr_sort.bam')


print ('Classifying reads between nuclear and mitochondrial sequences')
# Extract reads mapping to mitochondrial tRNAs.
os.system('samtools view -b '+sample_name+'_WGloc_only_trna_soft_clipped_rm_all_chr_sort.bam chrM > '+sample_name+'_WGloc_only_trna_mitochondrial.bam')
os.system('samtools sort '+sample_name+'_WGloc_only_trna_mitochondrial.bam'+ ' -o ' +sample_name+'_WGloc_only_trna_mitochondrial_sort.bam')
os.system('samtools index '+sample_name+'_WGloc_only_trna_mitochondrial_sort.bam')

# Remove reads mapping to mitochondrial tRNAs from the WG alignment files. 
classification_mitochondrial(sample_name)
os.system('picard FilterSamReads I='+sample_name+'_WGloc_only_trna_soft_clipped_rm_all_chr_sort.bam O='+sample_name+'_WGloc_only_trna_soft_clipped_rm.bam READ_LIST_FILE='+sample_name+'_mitochondrial_WG.txt FILTER=excludeReadList >/dev/null 2>&1')
os.system('samtools sort '+sample_name+'_WGloc_only_trna_soft_clipped_rm.bam'+ ' -o ' +sample_name+'_WGloc_only_trna_soft_clipped_rm_sort.bam')
os.system('samtools index '+sample_name+'_WGloc_only_trna_soft_clipped_rm_sort.bam')

print ('Classifying reads between precursor and mature sequences...')
# Extract the reads that contain leader,trailer and intronic regions (POSSIBLE PRECURSOR tRNAs).
os.system('samtools view -b -L ../../../Reference_Genomes/info/leader_trailer_intron.bed '+sample_name+'_WGloc_only_trna_soft_clipped_rm_sort.bam > '+sample_name+'_WGloc_only_trna_precursor_noFiltered.bam')
os.system('samtools sort '+sample_name+'_WGloc_only_trna_precursor_noFiltered.bam -o ' +sample_name+'_WGloc_only_trna_precursor_noFiltered_sort.bam')
os.system('samtools index '+sample_name+'_WGloc_only_trna_precursor_noFiltered_sort.bam')

# Filter reads that are mapping tRNAs that contain CCA sequences at the start of the trailing sequence.
classification_cca(sample_name, '_WGloc_only_trna_', tRNA_coordinates)

# Filter reads that contain a deletion in intronic regions.
classification_introns(sample_name,'_WGloc_only_trna_', tRNA_coordinates, intron_coordinates)
os.system('cat '+sample_name+'_exclude_precursor1.txt '+sample_name+'_exclude_precursor2.txt > '+sample_name+'_exclude_precursor.txt')

# Extract precursor reads.
with open(sample_name+'_exclude_precursor.txt') as file:
    first = file.read(1)
    
if not first:
    os.rename(sample_name+'_WGloc_only_trna_precursor_noFiltered_sort.bam',sample_name+'_WGloc_only_trna_precursor_Filtered.bam') 
else:    
    os.system('picard FilterSamReads I='+sample_name+'_WGloc_only_trna_precursor_noFiltered_sort.bam '+'O='+sample_name+'_WGloc_only_trna_precursor_Filtered.bam READ_LIST_FILE='+sample_name+'_exclude_precursor.txt FILTER=excludeReadList >/dev/null 2>&1')

os.system('samtools sort '+sample_name+'_WGloc_only_trna_precursor_Filtered.bam'+ ' -o ' +sample_name+'_WGloc_only_trna_precursor_Filtered_sort.bam')
os.system('samtools index '+sample_name+'_WGloc_only_trna_precursor_Filtered_sort.bam') 
os.system('samtools view -F 4 '+sample_name+'_WGloc_only_trna_precursor_Filtered_sort.bam '+'| cut -f1 | sort -u > '+sample_name+'_WGloc_only_trna_precursor_Filtered.txt')

# Extact mature reads.
os.system('picard FilterSamReads '+'I='+sample_name+'_WGloc_only_trna_soft_clipped_rm_sort.bam '+'O='+sample_name+'_WGloc_only_trna_mature.bam READ_LIST_FILE='+sample_name+'_WGloc_only_trna_precursor_Filtered.txt FILTER=excludeReadList >/dev/null 2>&1')
os.system('samtools sort '+sample_name+'_WGloc_only_trna_mature.bam'+ ' -o ' +sample_name+'_WGloc_only_trna_mature_sort.bam')
os.system('samtools index '+sample_name+'_WGloc_only_trna_mature_sort.bam')


# Merge the mature reads +  unmmaped reads, in order to peform the aligment with the mature genome:
os.system('samtools fastq '+sample_name+'_WGloc_only_trna_mature.bam > '+sample_name+'_WGloc_only_trna_mature.fastq')

# Merge fastq 
os.system('cat '+sample_name+'_WGloc_only_trna_mature.fastq '+sample_name+'_unmapped_WGloc.fastq > '+sample_name+'_WGloc_only_trna_mature_and_unmapped.fastq')

print ('\n')
