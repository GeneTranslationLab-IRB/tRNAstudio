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
import multiprocessing
import pysam
import pysamstats
import platform
import pandas as pd 
import Bio
from Bio import SeqIO
###############################################################################



def custom_removal_soft_clipped_bases(sample_name, file_name, tRNA_coordinates):
	''' Function used to perform a custom removal of soft clipped bases 
	tolerating 3’-CCA and 5’-G of tRNA_His'''

	# INPUT --> BAM file with soft clipped bases
	bamfile = pysam.AlignmentFile(sample_name+file_name+'sort.bam', "rb")
	
	# OUTPUT --> BAM file with soft clipped bases removed tolerating CCA sequences and and G addition to His
	out_bam = pysam.AlignmentFile(sample_name+file_name+'soft_clipped_rm_all_chr.bam', "wb", template=bamfile)

	tRNAsID = tRNA_coordinates['name'].tolist()

	readsIDs = {}

	for tRNAID in tRNAsID:	
		tRNAInfo = tRNA_coordinates.loc[tRNA_coordinates['name'] == tRNAID]

		if 'tRNAmt' in tRNAID:
			tRNAStart, tRNAEnd, tRNAChr   = int(list(tRNAInfo['Start'])[0]), int(list(tRNAInfo['End'])[0]), list(tRNAInfo['chrom'])[0]
		else:
			tRNAStart, tRNAEnd, tRNAChr   = int(list(tRNAInfo['Start'])[0]) + 50, int(list(tRNAInfo['End'])[0]) - 50, list(tRNAInfo['chrom'])[0]

		for read in bamfile.fetch(tRNAChr, int(list(tRNAInfo['Start'])[0]), int(list(tRNAInfo['End'])[0])):
			sequence = str(read.query_sequence)
			cigar_str = read.cigarstring
			if str(read.query_name) not in readsIDs:
				readsIDs[str(read.query_name)] = 1
				if 'S' in cigar_str:
					readStart, readEnd = read.reference_start, read.reference_end
					if read.is_reverse:
						if str(read.cigar[0])[1:-1].split(',')[0] == '4':
							soft_clipp=str(read.cigar[0])[1:-1].split(',')[1].replace(' ','')
							if int(soft_clipp) == 3 and sequence[:3] == 'TGG' and readStart == tRNAStart:
								pass
							elif int(soft_clipp) == 2 and sequence[:3] == 'TGG' and readStart+1 == tRNAStart:
								pass
							elif int(soft_clipp) == 1 and sequence[:3] == 'TGG' and readStart+2 == tRNAStart:
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
							elif int(soft_clipp) == 2 and sequence[-3:] == 'CCA' and readEnd-1 == tRNAEnd:
								pass
							elif int(soft_clipp) == 1 and sequence[-3:] == 'CCA' and readEnd-2 == tRNAEnd:
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



def custom_removal_soft_clipped_bases_2(sample_name,file_name):
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
	


def classification_cca(sample_name, file_name, tRNA_coordinates):
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
					


def classification_introns(sample_name, file_name, tRNA_coordinates, intron_coordinates):
	'''Function to classify reads that are mapping intronic regions'''

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


def classification_mitochondrial(sample):
	'''Function to filter mitochondrial tRNAs from the bam file containin both 
	processed and mitochondrial tRNAs'''

	bamfile_mitochondrial_reads_WG = pysam.AlignmentFile(sample+'_WGloc_only_trna_mitochondrial_sort.bam', "rb")
	mitochondrial_reads_WG = open(sample+'_mitochondrial_WG.txt',"w")
	for read in bamfile_mitochondrial_reads_WG.fetch():
		mitochondrial_reads_WG.write(str(read.query_name)+'\n')
	mitochondrial_reads_WG.close()

	# Mapping to a custom genome with only mitochondrial tRNA sequences.
	os.system('samtools fastq '+sample+'_WGloc_only_trna_mitochondrial_sort.bam > '+sample+'_WGloc_only_trna_mitochondrial.fastq')
	
	nucleos = multiprocessing.cpu_count()
	command = "bowtie2 --local -p "+ str(nucleos)
	os.system(command + ' -N 0 -x ../../../Reference_Genomes/mitochondrial_tRNA_refgenome '+sample+'_WGloc_only_trna_mitochondrial.fastq | samtools view -bSF4 - > ../Final_results/'+sample+'_mitochondrial.bam')
	os.system('samtools sort ../Final_results/'+sample+'_mitochondrial.bam -o ../Final_results/'+sample+'_mitochondrial_sort.bam')
	os.system('samtools index ../Final_results/'+sample+'_mitochondrial_sort.bam')

	bamfile_mitochondrial_reads_CG = pysam.AlignmentFile('../Final_results/'+sample+'_mitochondrial_sort.bam', "rb")
	mitochondrial_reads_CG = open(sample+'_mitochondrial_CG.txt',"w")

	for read in bamfile_mitochondrial_reads_CG.fetch():
		mitochondrial_reads_CG.write(str(read.query_name)+'\n')
	mitochondrial_reads_CG.close()

	os.system('rm *mitochondrial*.bam  *mitochondrial*.bai *mitochondrial*.fastq')

	os.system('picard FilterSamReads I='+sample+'_WGloc_only_trna_soft_clipped_rm_all_chr_sort.bam O='+sample+'_WGloc_mitochondrial_filtered.bam READ_LIST_FILE='+sample+'_mitochondrial_CG.txt FILTER=includeReadList >/dev/null 2>&1')
	os.system('samtools sort '+sample+'_WGloc_mitochondrial_filtered.bam -o '+sample+'_WGloc_mitochondrial_filtered_sort.bam')
	os.system('samtools index '+sample+'_WGloc_mitochondrial_filtered_sort.bam')




def counts(sample):
	'''Function used to obtain precursor, processed, mitochondrial 
	and the total number of counts  for each tRNA'''
	
	# Families IDs
	fam = open('../Reference_Genomes/info/tRNAsID.txt','r').readlines()
	fam = str(fam).replace('["','').replace('"]','')
	fam = ast.literal_eval(fam)

	sample = sample.split('.')[0]
	if not os.path.exists('../Results/'+sample+'/Counts/'):
		os.makedirs('../Results/'+sample+'/Counts/')
	os.chdir('../Results/'+sample+'/Counts/')

	# Annotation file 
	tRNAsAnnotation = '../../../Reference_Genomes/info/tRNAsAnnotation.gtf'
	
	# OUTPUTS --> Counts files
	os.system("featureCounts -t tRNAmat -a "+tRNAsAnnotation+' -o '+sample+"_counts_processed_tmp.txt ../Final_results/"+sample+"_processed_sort.bam")
	os.system("featureCounts -t tRNApre -a "+tRNAsAnnotation+' -o '+sample+"_counts_precursor_tmp.txt ../Final_results/"+sample+"_precursor_sort.bam")
	os.system("featureCounts -t tRNAmt -a "+tRNAsAnnotation+' -o '+sample+"_counts_mitochondrial_tmp.txt ../Final_results/"+sample+"_mitochondrial_sort.bam")

	# Read count files
	processed = pd.read_csv(sample+"_counts_processed_tmp.txt", comment='#', sep = '\t')
	precursor = pd.read_csv(sample+"_counts_precursor_tmp.txt", comment='#', sep = '\t')
	mitochondrial = pd.read_csv(sample+"_counts_mitochondrial_tmp.txt", comment='#', sep = '\t')   

	# Rename last column --> Counts
	precursor.columns = [*precursor.columns[:-1], 'Counts']
	processed.columns = [*processed.columns[:-1], 'Counts']
	mitochondrial.columns = [*mitochondrial.columns[:-1], 'Counts']

	precursorDef = precursor[['Geneid', 'Counts']]
	precursorDef.to_csv(sample+"_counts_precursor_tmp.txt", sep='\t', index = False, header = False)

	# Assign family to precursor tRNA
	for index, row in precursor.iterrows():
		tRNAid = row["Geneid"]
		tRNAfam = fam[tRNAid]
		precursor.loc[index, "Geneid"] = tRNAfam
	processedByFamCounts = processed[["Geneid", "Counts"]]
	# Group precursor counts by family
	precursorByFamCounts = precursor.groupby('Geneid')['Counts'].agg('sum')
	precursorByFamCounts = pd.DataFrame(precursorByFamCounts)
	precursorByFamCounts = precursorByFamCounts.reset_index()
	
	# Marge precursor and processed results
	total = pd.merge(precursorByFamCounts, processedByFamCounts, on="Geneid", how='inner')
	total['Counts'] = total['Counts_x'] + total['Counts_y']
	
	total,processed,mitochondrial  = total[['Geneid', 'Counts']], processed[['Geneid', 'Counts']], mitochondrial[['Geneid', 'Counts']]

	total.to_csv(sample+"_counts_total_tmp.txt", sep='\t', index = False, header = False)
	processed.to_csv(sample+"_counts_processed_tmp.txt", sep='\t', index = False, header = False )
	mitochondrial.to_csv(sample+"_counts_mitochondrial_tmp.txt", sep='\t', index = False, header = False )

	os.system('rm *summary*')



def mapping_quality(sample):
	'''Function used to obtain precursor, processed, mitochondrial 
	and the total number of counts  for each tRNA'''
	
	fam = open('../Reference_Genomes/info/tRNAsID.txt','r').readlines()
	fam = str(fam).replace('["','').replace('"]','')
	fam = ast.literal_eval(fam)
	
	os.chdir('../Results/'+sample+'/Counts/')

	# MAPQ analisys -- > mitochondrial 
	mitochondrialCounts = pd.read_csv(sample+"_counts_mitochondrial_tmp.txt", comment='#', sep = '\t', header = None)
	for index, row in mitochondrialCounts.iterrows():
		tRNA_ID = row[0]
		os.system('samtools view ../Final_results/'+sample+'_mitochondrial_sort.bam '+tRNA_ID+' | cut -f1 > '+tRNA_ID+'_IDs.txt')
		os.system('picard FilterSamReads I=../Alignment_WG/'+sample+'_WGloc_mitochondrial_filtered_sort.bam O='+sample+'_'+tRNA_ID+'.bam READ_LIST_FILE='+tRNA_ID+'_IDs.txt FILTER=includeReadList >/dev/null 2>&1')
		
		#Filter the reads with mapq  only include reads with mapping quality >= INT [0]. MAPQ >2 "GOOD MAPPING QUALITY READS"
		MitoReadsMappedMAPQ = int(subprocess.check_output('samtools view -c -bSq 3 '+sample+'_'+tRNA_ID+'.bam', shell=True).decode("utf-8").strip("\n"))
		os.system('rm *'+tRNA_ID+'*')
		mitochondrialCounts.loc[mitochondrialCounts[0] == tRNA_ID , 2] = str(MitoReadsMappedMAPQ)

	mitochondrialCounts['BadMAPQ'] = mitochondrialCounts[1].astype(int) - mitochondrialCounts[2].astype(int)
	mitochondrialCounts = mitochondrialCounts.drop(columns=[1])
	mitochondrialCounts.to_csv(sample+"_counts_mitochondrial.txt", sep='\t', index = False, header = False)

	# MAPQ analisys -- > Precursor 
	precursorCounts = pd.read_csv(sample+"_counts_precursor_tmp.txt", comment='#', sep = '\t', header = None)
	
	for index, row in precursorCounts.iterrows():
		tRNA_ID = row[0]
		PrecReadsMappedMAPQ = int(subprocess.check_output('samtools view -c -bSq 3 ../Final_results/'+sample+'_precursor_sort.bam '+tRNA_ID, shell=True).decode("utf-8").strip("\n"))
		precursorCounts.loc[precursorCounts[0] == tRNA_ID , 2] = str(PrecReadsMappedMAPQ)
	
	precursorCounts['BadMAPQ'] = precursorCounts[1].astype(int) - precursorCounts[2].astype(int)
	precursorCounts = precursorCounts.drop(columns=[1])
	precursorCounts.to_csv(sample+"_counts_precursor.txt", sep='\t', index = False, header = False)

	# MAPQ analisys -- > Processed 
	processedCounts = pd.read_csv(sample+"_counts_processed_tmp.txt", comment='#', sep = '\t', header = None)
	
	for index, row in processedCounts.iterrows():
		tRNA_ID = row[0]
		ProReadsMappedMAPQ = int(subprocess.check_output('samtools view -c -bSq 3 ../Final_results/'+sample+'_processed_sort.bam '+tRNA_ID, shell=True).decode("utf-8").strip("\n"))
		processedCounts.loc[processedCounts[0] == tRNA_ID , 2] = str(ProReadsMappedMAPQ)
	
	processedCounts['BadMAPQ'] = processedCounts[1].astype(int) - processedCounts[2].astype(int)
	processedCounts = processedCounts.drop(columns=[1])
	processedCounts.to_csv(sample+"_counts_processed.txt", sep='\t', index = False, header = False)

	# MAPQ analisys -- > Total
	# Assign family to precursor tRNA
	for index, row in precursorCounts.iterrows():
		tRNAid = row[0]
		tRNAfam = fam[tRNAid]
		precursorCounts.loc[index, 0] = tRNAfam
	
	# Group precursor counts by family
	precursorCounts[[2]] = precursorCounts[[2]].astype(int)
	precursorByFamCounts = precursorCounts.groupby(0)[2].agg('sum')
	precursorByFamCounts = pd.DataFrame(precursorByFamCounts)
	precursorByFamCounts = precursorByFamCounts.reset_index()
	processedByFamCounts = processedCounts[[0, 2]]


	totalCounts = pd.read_csv(sample+"_counts_total_tmp.txt", comment='#', sep = '\t', header = None)
	total = pd.merge(precursorByFamCounts, processedByFamCounts, on=0, how='inner')
	total.rename(columns={0: 'GeneID', '2_x': 'Counts_x', '2_y': 'Counts_y'}, inplace=True)
	total['Counts'] = total['Counts_x'].astype(int) + total['Counts_y'].astype(int)
	totalCounts['CountsMAPQ'] = total['Counts']
	totalCounts['BadMAPQ'] = totalCounts[1].astype(int) - totalCounts['CountsMAPQ'].astype(int)
	totalCounts = totalCounts.drop(columns=[1])
	totalCounts.to_csv(sample+"_counts_total.txt", sep='\t', index = False, header = False)
	os.system('rm *tmp*')
	os.system("cat "+sample+"_counts_mitochondrial.txt >> "+sample+"_counts_total.txt")




def pileup(sample):
	'''Pileup function to obtain base calls'''

	sample = sample.split('.')[0]

	if not os.path.exists('../Results/'+sample+'/Base_calling/'):
	    os.makedirs('../Results/'+sample+'/Base_calling/')

	os.chdir('../Results/'+sample+'/Base_calling')
	refgen_path = '../../../Reference_Genomes/'
	precursor_seq = SeqIO.parse(refgen_path+'precursor_fam_tRNA_refgenome.fa', "fasta")

	# Extract precursor length. 
	prec_length = {}
	for rec in precursor_seq:
		prec_length[rec.name] = len(rec.seq)

	# tRNA Coordinates.
	tRNA_coordinates = pd.read_csv(refgen_path+'info/tRNAsCoordinates.bed', sep='\t', header=None)
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
	fam=open(refgen_path+'info/tRNAsID.txt','r')

	# Create dictionary with nucleotide sequence for each family.
	family_sequence = SeqIO.to_dict(SeqIO.parse(refgen_path+'mature_fam_tRNA_refgenome.fa', "fasta"))
	fam=fam.readlines()
	fam=str(fam).replace('["','').replace('"]','')
	fam=ast.literal_eval(fam)

	# Correct positions reference (Additional bases)
	tRNA_pos_correct=open(refgen_path+'/info/tRNAsPositions.txt','r')
	tRNA_pos_ref={}
	for e in tRNA_pos_correct:
		e=e.split('\t')
		trna=e[0]
		tRNA_pos_ref[trna]=e[1].replace('\n','')

	# OUTPUT 
	tRNA_total=open(sample+'_base_calls_total.txt','w')
	tRNA_total.write('TRNA-POS'+'\t'+'A'+'\t'+'C'+'\t'+'G'+'\t'+'T'+'\t'+'REF-COUNTS'+'\n')

	# Input files (mature/processed and precursor)
	files=['../Final_results/'+sample+'_processed_sort.bam','../Final_results/'+sample+'_precursor_sort.bam']
	
	# Determine file type
	for file in files:

		bamfile = pysam.AlignmentFile(file)
		
		tRNA_all={}
		tRNA_all_prop={}
		
		bam_type = ''
		ref_genome = ''
		if 'processed' in file:
			bam_type='processed'
			ref_genome=refgen_path+'mature_fam_tRNA_refgenome.fa'
			file_base_call = open(sample+'_base_calls_processed.txt','w')
		if 'precursor' in file:
			bam_type='precursor'
			ref_genome = refgen_path+'precursor_fam_tRNA_refgenome.fa'
			file_base_call = open(sample+'_base_calls_precursor.txt','w')
		

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

			if bam_type=='processed':
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

						#cca_pos = ['74','75','76']
						#if not any(x in correct_pos for x in cca_pos):
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
	


