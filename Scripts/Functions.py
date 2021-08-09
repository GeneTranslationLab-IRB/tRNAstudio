#!/usr/bin/env python3

#### This script removes soft clipped bases without removing the CCA tail present in some reads as soft clipped bases ######

###############################################################################
# Built-in/Generic Imports and Modules #
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning) 
import os 
import re
import sys
import ast
import glob
import argparse
import subprocess
import pysam
import pysamstats
import platform
import pandas as pd 
import Bio
from Bio import SeqIO



###############################################################################


def CustomRemovalSoftClippedBases(sample_name, file_name, tRNA_coordinates):
	''' Function that removes the soft clipped bases tolerating the CCA tail and G addition to His'''

	# INPUT --> BAM file with soft clipped bases
	bamfile = pysam.AlignmentFile(sample_name+file_name+'sort.bam', "rb")
	
	# OUTPUT --> BAM file with soft clipped bases removed tolerating CCA sequences and and G addition to His
	out_bam = pysam.AlignmentFile(sample_name+file_name+'soft_clipped_rm.bam', "wb", template=bamfile)

	precursor_seq = SeqIO.parse('../../../Reference_Genomes/precursor_tRNA_refgenome.fa', "fasta")

	for rec in precursor_seq:	
		tRNAID = rec.name
		tRNAInfo = tRNA_coordinates.loc[tRNA_coordinates['name'] == (tRNAID.replace('(+)',('')).replace('(-)',('')))]
		tRNAStart, tRNAEnd, tRNAChr   = int(list(tRNAInfo['Start'])[0]) + 50, int(list(tRNAInfo['End'])[0]) - 50, list(tRNAInfo['chrom'])[0]

		for read in bamfile.fetch(tRNAChr, int(list(tRNAInfo['Start'])[0]), int(list(tRNAInfo['End'])[0])):
			sequence = str(read.query_sequence)
			cigar_str = read.cigarstring

			if 'S' in cigar_str:
				readStart, readEnd = read.reference_start, read.reference_end

				if read.is_reverse:
					if str(read.cigar[0])[1:-1].split(',')[0] == '4':
						soft_clipp=str(read.cigar[0])[1:-1].split(',')[1].replace(' ','')
						if int(soft_clipp) == 3 and sequence[:3] == 'TGG' and readStart == tRNAStart:
							pass
						elif int(soft_clipp) == 2 and sequence[:2] == 'GG' and readStart == tRNAStart:
							pass
						elif int(soft_clipp) == 1 and sequence[:1] == 'G' and readStart == tRNAStart:
							pass
						elif int(soft_clipp) == 2 and sequence[:3] == 'TGG' and readStart+1 == tRNAStart:
							pass
						elif int(soft_clipp) == 1 and sequence[:3] == 'TGG' and readStart+2 == tRNAStart:
							pass
						elif int(soft_clipp) == 1 and sequence[:2] == 'GG' and readStart+1 == tRNAStart:
							pass

						else:
							soft_clipp_seq=str(sequence[:int(soft_clipp)])
							if soft_clipp_seq[-3:] == 'TGG' and int(soft_clipp) > 3 and readStart == tRNAStart:
								
								read.cigarstring = '3S'+str(cigar_str[(len(soft_clipp)+1):])
								qual = read.query_qualities[(int(soft_clipp)-3):]
								read.query_sequence , read.query_qualities = str(sequence[(int(soft_clipp)-3):]), qual
								cigar_str, sequence = read.cigarstring, str(read.query_sequence)

							else:
								qual = read.query_qualities[int(soft_clipp):]
								read.cigarstring = cigar_str[(len(soft_clipp)+1):]
								read.query_sequence, read.query_qualities = sequence[int(soft_clipp):], qual
								cigar_str, sequence = read.cigarstring, str(read.query_sequence)

					if str(read.cigar[-1])[1:-1].split(',')[0] == '4':
						soft_clipp = str(read.cigar[-1])[1:-1].split(',')[1].replace(' ','')
						soft_clipp_seq=str(sequence[-int(soft_clipp):])

						if 'His' in tRNAID and soft_clipp_seq == 'C' and readEnd == tRNAEnd and int(soft_clipp) == 1:
							out_bam.write(read)
						
						elif 'His' in tRNAID and soft_clipp_seq[0] == 'C' and readEnd == tRNAEnd and int(soft_clipp) > 1:
							read.cigarstring = str(cigar_str[:-(len(soft_clipp)+1)])+'1S'
							qual = read.query_qualities[:-(int(soft_clipp)-1)]
							read.query_sequence, read.query_qualities = str(sequence[:-(int(soft_clipp)-1)]), qual
							cigar_str, sequence = read.cigarstring, str(read.query_sequence)
							out_bam.write(read)

						else:
							read.cigarstring = cigar_str[:-(len(soft_clipp)+1)]
							qual = read.query_qualities[:-int(soft_clipp)]
							read.query_sequence, read.query_qualities = str(sequence[:-int(soft_clipp)]), qual
							cigar_str = read.cigarstring
							out_bam.write(read)

					else:
						out_bam.write(read)

				else:				
					if str(read.cigar[-1])[1:-1].split(',')[0] == '4':
						soft_clipp=str(read.cigar[-1])[1:-1].split(',')[1].replace(' ','')
						if int(soft_clipp) == 3 and sequence[-3:] == 'CCA' and readEnd == tRNAEnd:
							pass
						elif int(soft_clipp) == 2 and sequence[-2:] == 'CC' and readEnd == tRNAEnd:
							pass
						elif int(soft_clipp) == 1 and sequence[-1:] == 'C' and readEnd == tRNAEnd:
							pass
						elif int(soft_clipp) == 2 and sequence[-3:] == 'CCA' and readEnd-1 == tRNAEnd:
							pass
						elif int(soft_clipp) == 1 and sequence[-3:] == 'CCA' and readEnd-2 == tRNAEnd:
							pass
						elif int(soft_clipp) == 1 and sequence[-2:] == 'CC' and readEnd-1 == tRNAEnd:
							pass

						else:
							soft_clipp_seq=str(sequence[-int(soft_clipp):])
							if soft_clipp_seq[:3] == 'CCA' and int(soft_clipp) > 3 and  readEnd == tRNAEnd:
								read.cigarstring = str(cigar_str[:-(len(soft_clipp)+1)])+'3S'
								qual = read.query_qualities[:-(int(soft_clipp)-3)]
								read.query_sequence, read.query_qualities = str(sequence[:-(int(soft_clipp)-3)]), qual
								cigar_str, sequence = read.cigarstring, str(read.query_sequence)

							else:
								read.cigarstring = cigar_str[:-(len(soft_clipp)+1)]
								qual = read.query_qualities[:-int(soft_clipp)]
								read.query_sequence, read.query_qualities = str(sequence[:-int(soft_clipp)]), qual
								cigar_str, sequence = read.cigarstring, str(read.query_sequence)

					if str(read.cigar[0])[1:-1].split(',')[0] == '4':
						soft_clipp = str(read.cigar[0])[1:-1].split(',')[1].replace(' ','')
						soft_clipp_seq=str(sequence[:int(soft_clipp)])
						
						if 'His' in tRNAID and soft_clipp_seq == 'G' and readStart == tRNAStart and int(soft_clipp) == 1:
							out_bam.write(read)

						elif 'His' in tRNAID and soft_clipp_seq[-1:] == 'G' and readStart == tRNAStart and int(soft_clipp) > 1:
							read.cigarstring = '1S'+str(cigar_str[(len(soft_clipp)+1):])
							qual = read.query_qualities[(int(soft_clipp)-1):]
							read.query_sequence , read.query_qualities = str(sequence[(int(soft_clipp)-1):]), qual
							cigar_str, sequence = read.cigarstring, str(read.query_sequence)
							out_bam.write(read)

						else:
							qual = read.query_qualities[int(soft_clipp):]
							read.cigarstring = cigar_str[(len(soft_clipp)+1):]
							read.query_sequence, read.query_qualities = sequence[int(soft_clipp):], qual
							cigar_str, sequence = read.cigarstring, str(read.query_sequence)
							out_bam.write(read)
							
					else:
						out_bam.write(read)
			else:
				out_bam.write(read)

	out_bam.close()
	bamfile.close()




def CustomRemovalSoftClippedBases2(sample_name,file_name):
	''' Remove soft clipped bases '''
	
	# INPUT --> file with soft clipped bases.
	bamfile = pysam.AlignmentFile(sample_name+file_name+'sort.bam', "rb")
	
	# OUTPUT --> file with soft clipped bases removed.
	out_bam = pysam.AlignmentFile(sample_name+file_name+'soft_clipped_rm.bam', "wb", template=bamfile)
	
	for read in bamfile.fetch():
		cigar_str = read.cigarstring
		
		if 'S' in cigar_str:
			sequence = str(read.query_sequence)
			soft_clipp = str(read.cigar[0])[1:-1].split(',')[1].replace(' ','')
			if str(read.cigar[0])[1:-1].split(',')[0] == '4':
				qual = read.query_qualities[int(soft_clipp):]
				read.cigarstring = cigar_str[(len(soft_clipp)+1):]
				read.query_sequence, read.query_qualities = sequence[int(soft_clipp):], qual
				cigar_str, sequence = read.cigarstring, str(read.query_sequence)
			
			if str(read.cigar[-1])[1:-1].split(',')[0] == '4':
				soft_clipp = str(read.cigar[-1])[1:-1].split(',')[1].replace(' ','')
				qual = read.query_qualities[:-int(soft_clipp)]
				read.cigarstring = cigar_str[:-(len(soft_clipp)+1)]
				read.query_sequence, read.query_qualities = str(sequence[:-int(soft_clipp)]), qual
				cigar_str, sequence = read.cigarstring, str(read.query_sequence)

				out_bam.write(read)	

			else:
				out_bam.write(read)	
		else:
			out_bam.write(read)	

	out_bam.close()
	bamfile.close()
	



def ClassificationCCA(sample_name, file_name, tRNA_coordinates):
	'''Function that will take into account reads that are mapping tRNAs 
	that contain CCA sequences at the start of the trailer sequence'''

	# INPUT --> BAM file with "precursor" sequences to filter.
	bamfilePrecursorFilter = pysam.AlignmentFile(sample_name+file_name+'precursor_noFiltered_sort.bam', "rb")

	# OUTPUT --> Reads classified as mature/"processed".
	mature_reads = open(sample_name+'_exclude_precursor2.txt',"w")
	precursor_seq = SeqIO.parse('../../../Reference_Genomes/precursor_tRNA_refgenome.fa', "fasta")

	for rec in precursor_seq:
		trailer=str(rec.seq)[-50:].upper()
		tRNAID=rec.name
		tRNAInfo = tRNA_coordinates.loc[tRNA_coordinates['name'] == (tRNAID.replace('(+)',('')).replace('(-)',('')))]
		tRNAStart, tRNAEnd, tRNAChr = int(list(tRNAInfo['Start'])[0])+50, int(list(tRNAInfo['End'])[0])-50, list(tRNAInfo['chrom'])[0]
		
		if trailer.startswith('C'):
			for read in bamfilePrecursorFilter.fetch(tRNAChr, tRNAStart, tRNAEnd):
				cigar_str=read.cigarstring
				if 'S' not in cigar_str:
					readStart, readEnd = (read.reference_start), (read.reference_end)
					possible_cca_seq=read.query_sequence[-(readEnd-tRNAEnd):]

					if not read.is_reverse and possible_cca_seq == 'CCA' or possible_cca_seq == 'CC' or possible_cca_seq == 'C':
						if 1 <= readEnd-tRNAEnd <= 3 and readStart>=tRNAStart:						
							mature_reads.write(read.query_name+'\n')

					possible_cca_seq_rev=read.query_sequence[:tRNAStart-readStart]	
					if read.is_reverse and possible_cca_seq_rev  == 'TGG' or possible_cca_seq_rev == 'GG' or possible_cca_seq_rev == 'G':
						if 1 <= tRNAStart-readStart <= 3 and readEnd<=tRNAEnd:
							mature_reads.write(read.query_name+'\n')

		if trailer.startswith('GCA') or trailer.startswith('TCA') or trailer.startswith('ACA'):
			for read in bamfilePrecursorFilter.fetch(tRNAChr, tRNAStart, tRNAEnd):
				cigar_str=read.cigarstring
				
				if 'S' not in cigar_str:
					readStart, readEnd = read.reference_start, read.reference_end
					possible_cca_seq=read.query_sequence[-(readEnd-tRNAEnd):]
					
					if not read.is_reverse and possible_cca_seq == 'CCA':
						if readEnd-tRNAEnd == 3:				
							mature_reads.write(read.query_name+'\n')
	
					possible_cca_seq_rev=read.query_sequence[:tRNAStart-readStart]	
					if read.is_reverse and possible_cca_seq_rev  == 'TGG':
						if tRNAStart-readStart == 3:
							mature_reads.write(read.query_name+'\n')

	for read in bamfilePrecursorFilter.fetch():
		cigar_str=read.cigarstring
		if 'S' in cigar_str:
			mature_reads.write(read.query_name+'\n')
	
	mature_reads.close()
					



def ClassificationIntrons(sample_name, file_name, tRNA_coordinates, intron_coordinates):
	''' Function to classify reads that are mapping intronic regions '''

	# INPUT --> BAM file with "precursor" sequences to filter.
	bamfilePrecursorFilter = pysam.AlignmentFile(sample_name+file_name+'precursor_noFiltered_sort.bam', "rb")
	
	# OUTPUT --> Reads classified as mature/"processed".
	mature_reads = open(sample_name+'_exclude_precursor1.txt',"w")
	
	for index, row in intron_coordinates.iterrows():

		IntronStart, IntronEnd, tRNAChr, tRNAID= int(row['Start']), int(row['End']), row['chrom'], row['name']
		IntronLen = IntronEnd-IntronStart

		tRNAInfo = tRNA_coordinates.loc[tRNA_coordinates['name'] == tRNAID]
		tRNAStart, tRNAEnd = int(list(tRNAInfo['Start'])[0])+50, int(list(tRNAInfo['End'])[0])-50

		for read in bamfilePrecursorFilter.fetch(tRNAChr, IntronStart , IntronEnd):
			for e in read.cigartuples:
				if 'S' in read.cigarstring:
					mature_reads.write(read.query_name+'\n')
				else:
					d=0
					for e in read.cigartuples:
						if int(list(e)[0]) == 2 and int(list(e)[1]) == IntronLen:
							d=+1
					if d > 0:
						mature_reads.write(read.query_name+'\n')





def counts(sample):
	
	sample = sample.split('.')[0]
	if not os.path.exists('../Results/'+sample+'/Counts/'):
		os.makedirs('../Results/'+sample+'/Counts/')
	
	os.chdir('../Results/'+sample+'/Counts/')

	print ('Obtaining counts ...')

	# Mature
	os.system('samtools idxstats '+'../Final_results/'+sample+'all_mature_sort.bam >'+sample+'all_mature_sort.txt')
	mature_counts=open(sample+'all_mature_sort.txt','r')
	
	# Precursor
	os.system('samtools idxstats '+'../Final_results/'+sample+'_PGloc_mapped_sort.bam >'+sample+'all_precursor_sort.txt')
	precursor_counts=open(sample+'all_precursor_sort.txt','r')

	# tRNA family ID
	fam=open('../../../Reference_Genomes/info/families_grups/IDs_tRNA.txt','r')
	fam=fam.readlines()
	fam=str(fam).replace('["','').replace('"]','')
	fam=ast.literal_eval(fam)


	total_counts={}

	for trna in mature_counts:
		trna=trna.split('\t')
		total_counts[trna[0]]=trna[2]

	for trna in precursor_counts:
		trna=trna.split('\t')
		trna_id=str(trna[0])
		if trna_id.startswith('tRNA'):
			fam_trna=fam[trna_id]
			if fam_trna in total_counts:
				total_counts[fam_trna]=int(total_counts[fam_trna])+int(trna[2])

	# Total
	total=open(sample+'all_total.txt','w')
	for trna in sorted(total_counts):
		total.write(trna+'\t'+str(total_counts[trna])+'\n')





def pileup(sample):

	sample = sample.split('.')[0]

	if not os.path.exists('../Results/'+sample+'/Base_calling/'):
	    os.makedirs('../Results/'+sample+'/Base_calling/')

	os.chdir('../Results/'+sample+'/Base_calling')
	

	refgen_path = '../../../Reference_Genomes/'
	
	precursor_seq = SeqIO.parse(refgen_path+'precursor_fam_tRNA_refgenome.fa', "fasta")

	print ('Obtaining pileup format ...')


	# Extract precursor length. 
	prec_length = {}
	for rec in precursor_seq:
		prec_length[rec.name] = len(rec.seq)


	# tRNA Coordinates.
	tRNA_coordinates = pd.read_csv(refgen_path+'info/tRNA_only.bed', sep='\t', header=None)
	header = ['chrom', 'Start', 'End', 'name', 'score', 'strand']
	tRNA_coordinates.columns = header[:len(tRNA_coordinates.columns)]

	# Extract intron info.
	intron_coordinates = pd.read_csv(refgen_path+'info/intron.bed', sep='\t', header=None)
	header = ['chrom', 'Start', 'End', 'name', 'strand']
	intron_coordinates.columns = header[:len(intron_coordinates.columns)]

	tRNA_intron = {}

	for index, row in intron_coordinates.iterrows():
		IntronStart, IntronEnd, tRNAID_1, strand= int(row['Start']), int(row['End']), row['name'], row['strand']
		for index, row in tRNA_coordinates.iterrows():
			tRNAStart, tRNAEnd, tRNAID_2 = int(row['Start'])+50, int(row['End']-50), row['name']

			if tRNAID_1 == tRNAID_2:
				if strand == '+':
					Start = IntronStart-tRNAStart
					End = (IntronEnd-IntronStart) + (IntronStart-tRNAStart)
					Length = IntronEnd-IntronStart
					new_row = str(Start)+' '+str(End)+' '+str(Length)
					tRNA_intron[tRNAID_1] = new_row.split(' ')

				
				if strand == '-':
					Start =  tRNAEnd-IntronEnd
					End = ( tRNAEnd-IntronEnd) + (IntronEnd-IntronStart)
					Length = IntronEnd-IntronStart
					new_row = str(Start)+' '+str(End)+' '+str(Length)
					tRNA_intron[tRNAID_1] = new_row.split(' ')





	# Dictionary with tRNAs IDs 
	fam=open(refgen_path+'info/families_grups/IDs_tRNA.txt','r')

	# Create dictionary with nucleotide sequence for each family.
	family_sequence = SeqIO.to_dict(SeqIO.parse(refgen_path+'families_tRNA_refgenome.fa', "fasta"))
	fam=fam.readlines()
	fam=str(fam).replace('["','').replace('"]','')
	fam=ast.literal_eval(fam)


	# Correct positions reference (Additional bases)
	tRNA_pos_correct=open(refgen_path+'/Modifications_files/Ref_seq/trna_pos.txt','r')
	tRNA_pos_ref={}
	for e in tRNA_pos_correct:
		e=e.split('\t')
		trna=e[0]
		tRNA_pos_ref[trna]=e[1].replace('\n','')


	# OUTPUT 
	tRNA_total=open(sample+'_all_base_calling_by_pos_CORRECT_OK.txt','w')
	tRNA_total.write('TRNA-POS'+'\t'+'A'+'\t'+'C'+'\t'+'G'+'\t'+'T'+'\t'+'REF-COUNTS'+'\n')


	# Input files (mature/processed and precursor)
	files=['../Final_results/'+sample+'all_mature_sort.bam','../Final_results/'+sample+'_PGloc_mapped_sort.bam']


	tRNA_intron_info=open(refgen_path+'info/intron_info.txt','r')



	# Determine file type
	for file in files:

		bamfile = pysam.AlignmentFile(file)
		
		tRNA_all={}
		tRNA_all_prop={}
		
		bam_type=''
		ref_genome=''
		if 'mature' in file:
			bam_type='mature'
			ref_genome=refgen_path+'families_tRNA_refgenome.fa'
		if 'PG' in file:
			bam_type='precursor'
			ref_genome=refgen_path+'precursor_fam_tRNA_refgenome.fa'
		
		file_base_call=open(sample+'_'+bam_type+'_'+'base_calling_CORRECT_OK.txt','w')
		

		# BASE CALLING
		for record in pysamstats.stat_variation(bamfile,  fafile=ref_genome,  max_depth=10000000):
			
			ref_base = record['ref']
			tRNA_info = 'REF-'+str(record['ref'])+':'+str(record[ref_base])+' '+'A:'+str(record['A'])+' '+'C:'+str(record['C'])+' '+'G:'+str(record['G'])+' '+'T:'+str(record['T'])
			tRNA = record['chrom']
			pos = ''
			ref = record['ref']

			if bam_type=='precursor':
				pos_i=record['pos']
				if pos_i > 49 and pos_i < int(prec_length[tRNA])-50:
					pos = pos_i - 49
					if tRNA in tRNA_intron:
						info_intron=tRNA_intron[tRNA]
						start_intron=int(info_intron[0])
						end_intron=int(info_intron[1])
						length_intron=int(info_intron[2])

						if pos > start_intron and pos <= end_intron:
							pos=''
						else:
							if pos > end_intron:
								pos=pos-length_intron
				
					if 'His' in tRNA:
						pos=pos+1

			if bam_type=='mature':
				pos=record['pos']+1

			tRNA_info=tRNA_info.split(' ')
			tRNA_info.remove(str(ref_base)+':'+str(record[ref_base]))
			if pos !='':
				if tRNA in tRNA_all:
					tRNA_all[tRNA][pos]=tRNA_info

				else:
					tRNA_all[tRNA]={}
					tRNA_all[tRNA][pos]=tRNA_info
			
			total_bases=0
				
			
		for trna in tRNA_all:
			file_base_call.write('>'+trna+'\n')
			positions=tRNA_all[trna]
			for position in positions:
				file_base_call.write(str(position)+'\t'+str(positions[position])+'\n')

		for trna in tRNA_all:
			
			positions=tRNA_all[trna]
			pos_used=[]

			for position in positions:
				bases=positions[position]
				def take_first(elem):
					return elem.split(':')[0][-1]
				
				sorted_bases = sorted(bases, key = take_first)
				sorted_bases=','.join(sorted_bases)
				sorted_bases=sorted_bases.split(',')
				values=''
				ref=''
				for e in bases:
					if 'REF' in e:
						ref=e
				for val in sorted_bases:
					val=val.split(':')[1]
					values=values+val+'\t'
				values=values[:-1]
				
				if bam_type=='precursor':
					trna_fam=fam[trna]
					if trna_fam in tRNA_pos_ref:
						correct_pos=tRNA_pos_ref[trna_fam]

						#In the dictionary tRNA pos ref we have for each trna the positions of reference. But we have to take in to account the addition of  the cca tail. That is why we have to add the positions 74,75,76.
						pos=(correct_pos+','+'74'+','+'75'+','+'76').replace('\r','')
						pos=pos.split(',')
						correct_pos=pos[position-1]
						pos_used.append(correct_pos)
						tRNA_total.write(trna_fam+':'+str(position)+':'+ref.split(':')[0]+':'+correct_pos+'\t'+values+'\t'+ref.split(':')[1]+'\n')

				else:
					if trna in tRNA_pos_ref:
						correct_pos=tRNA_pos_ref[trna]
						

						pos=(correct_pos+','+'74'+','+'75'+','+'76').replace('\r','')
						pos=pos.split(',')
						correct_pos=pos[position-1]
						pos_used.append(correct_pos)
						tRNA_total.write(trna+':'+str(position)+':'+ref.split(':')[0]+':'+correct_pos+'\t'+values+'\t'+ref.split(':')[1]+'\n')
			

			#Positions with value 0
		
			if bam_type=='precursor':
				trna=fam[trna]
				pos_ref = tRNA_pos_ref[trna].replace('\r','')
				pos_ref = pos_ref.split(',')
				pos = list(range(1, len(pos_ref)+1))
				for e in range(len(pos_ref)):
					if pos_ref[e] not in pos_used:
						values='0'+'\t'+'0'+'\t'+'0'+'\t'+'0'
						ref_nuc_seq=(family_sequence[trna].seq)[e].upper()
						tRNA_total.write(trna+':'+str(pos[e])+':REF-'+ref_nuc_seq+':'+pos_ref[e]+'\t'+values+'\t'+'0'+'\n')
				

			else:
				pos_ref = tRNA_pos_ref[trna].replace('\r','')
				pos_ref = pos_ref.split(',')			
				pos = list(range(1, len(pos_ref)+1))
				for e in range(len(pos_ref)):
					if pos_ref[e] not in pos_used:
						values='0'+'\t'+'0'+'\t'+'0'+'\t'+'0'
						ref_nuc_seq=(family_sequence[trna].seq)[e].upper()
						tRNA_total.write(trna+':'+str(pos[e])+':REF-'+ref_nuc_seq+':'+pos_ref[e]+'\t'+values+'\t'+'0'+'\n')
	
		

	