
import sys
import os
import subprocess
#Mapping 

sample_name = 'SRR1836127'
os.chdir('../Results/'+sample_name)

#os.system('samtools view -c Alignment_WG/'+sample_name+'_WGloc_mapped.bam')


readsMapped = subprocess.check_output('samtools view -c Alignment_WG/'+sample_name+'_WGloc_mapped_sort.bam', shell=True)

readstRNA = subprocess.check_output('samtools view -c Alignment_WG/'+sample_name+'_WGloc_only_trna_sort.bam', shell=True)

readsPrecursor = subprocess.check_output('samtools view -c Alignment_WG/'+sample_name+'_WGloc_only_trna_precursor_Filtered_sort.bam', shell=True)

readsMature = subprocess.check_output('samtools view -c Alignment_WG/'+sample_name+'_WGloc_only_trna_mature_sort.bam', shell=True)

readsMatureMG = subprocess.check_output('samtools view -c Alignment_MG/'+sample_name+'_MGloc_mapped_sort.bam', shell=True)

readsPrecursorPG = subprocess.check_output('samtools view -c Alignment_PG/'+sample_name+'_PGloc_mapped_sort.bam', shell=True)

readsMatureM1G = subprocess.check_output('samtools view -c Alignment_M1G/'+sample_name+'_MG_1M_loc_mapped_sort.bam', shell=True)
'''

test = subprocess.check_output('samtools view -c -q 3 Final_results/'+sample_name+'all_mature_sort.bam', shell=True)

print (test)

test = subprocess.check_output('samtools view -c -q 3 Final_results/'+sample_name+'_PGloc_mapped_sort.bam', shell=True)

print (test)
'''
print (sample_name)
print ('reads mapped WG')
print (readsMapped)
print ('reads Mapped tRNAs WG')
print (readstRNA)
print ('Reads Precursor WG')
print (readsPrecursor)
print ('Reads Mature WG')
print (readsMature)
print ('Reads Mature MG')
print (readsMatureMG)
print ('reads Mature M1G')
print (readsMatureM1G)
print ('Reads Precursor PG')
print (readsPrecursorPG)
