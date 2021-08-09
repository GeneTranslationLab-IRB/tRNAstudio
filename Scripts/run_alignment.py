#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 13:28:22 2020

@author: ignacio
"""
import os
import subprocess
import sys 

def run_alignment(file, folder):
    '''
    This function calls the pipeline :)
    '''
    scripts_folder = folder[:-17] + "/Scripts"
    os.chdir(scripts_folder)
    print("Starting the pipeline.")
    print("Checking for missing modules and installing them if missing.")
    os.system('python3 modules.py')
    print("Done.")
    print("Starting the whole genome alignment. This takes time.")
    
    os.system('python3 Aln_WG.py '+file+' '+folder)
    print("Done!")
    print("Aligning versus the mature genome.")
    os.system('python3 Aln_MG.py '+file)
    print("Done!")
    print("Aligning versus the precursor genome.")
    os.system('python3 Aln_PG.py '+file)
    print("Done!")
    print("Aligning versus the mature genome with one mismatch in the seed.")
    os.system('python3 Aln_M1G.py '+file)
    print("Done!")
    print("Obtaining the final counts.")
    os.system('python3 Obtain_counts.py '+file)
    print("Done!")
    print("Doing the pileup.")
    os.system('python3 pileup_mod.py '+file)
    print("Done!")
    
    
cwd = os.getcwd()


fastq_directory = cwd[:-8] + "/Fastq_downloaded/"
"""
os.chdir(fastq_directory)

file = open("sample_data.txt", "r")
sample_data = file.readlines()
sample_data = sample_data[1:]
for sample in sample_data:
    sample = sample.split("\t")[0]
    sample = sample + ".fastq"
    run_alignment(sample, fastq_directory)
"""
sample = sys.argv[1] + ".fastq"
run_alignment(sample, fastq_directory)
    

print("Runing R preparation")
#subprocess.call(["/usr/bin/Rscript", "--vanilla", "Join_results.r"])
print("Done")

    
    
    
