#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 09:10:44 2020

@author: ignacio
"""

import subprocess
import os

def r_preparation():
    '''
    This script takes all the counts step result files and joins them in a single file in order to analyse them in R
    '''
    os.chdir("../Results/R_files")
    if "Base_calling" not in os.listdir(os.getcwd()):    
        os.mkdir("Base_calling")
    if "Counts" not in os.listdir(os.getcwd()):
        os.mkdir("Counts")
    if "Mapping_quality" not in os.listdir(os.getcwd()):
        os.mkdir("Mapping_quality")
        
    file = open("../../Fastq_downloaded/sample_data.txt", "r")
    samples = []
    for line in file.readlines():
        samples.append(line.split("\t")[0])
    file.close()
    samples = samples[1:]  # Remove the header.
        
    if "Linux" in o_sys or "Darwin" in o_sys:
        os.chdir("..")
        for sample in samples:
            path_to_file1 = os.getcwd()+"/"+sample+"/Counts/"+sample+\
                "_counts_total.txt" 
            path_to_file2 = os.getcwd()+"/"+sample+"/Counts/"+sample+\
                "_counts_processed.txt"                 
            path_to_file3 = os.getcwd()+"/"+sample+"/Base_calling/"+sample+ \
                "_base_calls_total.txt" 
            path_to_file4 = os.getcwd()+"/"+sample+"/Alignment_WG/"+sample+\
                "_WGloc_only_trna.bam"
            path_to_file5 = os.getcwd()+"/"+sample+"/Alignment_MG/"+sample+\
                "_MGloc_mapped.bam"
            path_to_file6 = os.getcwd()+"/"+sample+"/Alignment_PG/"+sample+\
                "_precursor.bam"
            path_to_file7 = os.getcwd()+"/"+sample+"/Alignment_M1G/"+sample+\
                "_MG_1M_loc_mapped.bam"            

            cmd1 = "cp " + path_to_file1 + " " + os.getcwd() + "/R_files" + \
                "/Counts"
            cmd2 = "cp " + path_to_file2 + " " + os.getcwd() + "/R_files" + \
                "/Counts"
            cmd3 = "cp " + path_to_file3 + " " + os.getcwd() + "/R_files" + \
                "/Base_calling"
            cmd4 = "cp " + path_to_file4 + " " + os.getcwd() + "/R_files" + \
                "/Mapping_quality"
            cmd5 = "cp " + path_to_file5 + " " + os.getcwd() + "/R_files" + \
                "/Mapping_quality"
            cmd6 = "cp " + path_to_file6 + " " + os.getcwd() + "/R_files" + \
                "/Mapping_quality"
            cmd7 = "cp " + path_to_file7 + " " + os.getcwd() + "/R_files" + \
                "/Mapping_quality"
                
            os.system(cmd1)
            os.system(cmd2)
            os.system(cmd3)
            os.system(cmd4)
            os.system(cmd5)
            os.system(cmd6)
            os.system(cmd7)
        os.chdir("../Scripts")
        subprocess.call(["/usr/bin/Rscript", "--vanilla", "Join_results.r"])
        
    if "Windows" in o_sys:  # This is gonna fail.
        os.chdir("..")
        for sample in samples:
            path_to_file1 = os.getcwd()+"\ "+sample+"\Counts\ "+sample+\
                "_counts_total.txt" 
            path_to_file2 = os.getcwd()+"\ "+sample+"\Counts\ "+sample+\
                "_counts_processed.txt"    
            path_to_file3 = os.getcwd()+"\ "+sample+"\Base_calling\ "+sample+ \
                "_base_calls_total.txt" 
            path_to_file4 = os.getcwd()+"\ "+sample+"\Alignment_WG\ "+sample+\
                "_WGloc_only_trna.bam"
            path_to_file5 = os.getcwd()+"\ "+sample+"\Alignment_MG\ "+sample+\
                "_MGloc_mapped.bam"
            path_to_file6 = os.getcwd()+"\ "+sample+"\Alignment_PG\ "+sample+\
                "_precursor.bam"
            path_to_file7 = os.getcwd()+"\ "+sample+"/Alignment_M1G\ "+sample+\
                "_MG_1M_loc_mapped.bam"  
                
            cmd1 = "copy " + path_to_file1 + " " + os.getcwd() + "\R_files" + \
                "\Counts"
            cmd2 = "copy " + path_to_file1 + " " + os.getcwd() + "\R_files" + \
                "\Counts"
            cmd3 = "copy " + path_to_file3 + " " + os.getcwd() + "\R_files" + \
                "\Base_calling"
            cmd4 = "cp " + path_to_file4 + " " + os.getcwd() + "\R_files" + \
                "\Mapping_quality"
            cmd5 = "cp " + path_to_file5 + " " + os.getcwd() + "\R_files" + \
                "\Mapping_quality"
            cmd6 = "cp " + path_to_file6 + " " + os.getcwd() + "\R_files" + \
                "\Mapping_quality"
            cmd7 = "cp " + path_to_file7 + " " + os.getcwd() + "\R_files" + \
                "\Mapping_quality"
                
            os.system(cmd1)
            os.system(cmd2)    
            os.system(cmd3)
            os.system(cmd4)
            os.system(cmd5)
            os.system(cmd6)
            os.system(cmd7)
        
        
o_sys = "Linux"
r_preparation()