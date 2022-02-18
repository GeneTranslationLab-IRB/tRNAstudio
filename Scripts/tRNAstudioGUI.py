#!/usr/bin/env python3

'''
Pipeline for the analisis of tRNA-Seq Datasets (tRNAstudio GUI)
'''


###############################################################################
print ('\n'+'############# tRNAstudio ###############'+'\n')

# Built-in/Generic Imports #
import sys
import os
import platform
import subprocess
from subprocess import check_output
import pysam
import tkinter 
from tkinter import messagebox
from tkinter import filedialog
from tkinter import ttk
from tkinter import font
import tkinter.font as tkFont
from tkinter import messagebox as mb
from time import time
import multiprocessing
import pandas as pd
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

# Import functions #
from Functions import counts
from Functions import pileup
from Functions import mapping_quality

# System
o_sys=platform.system()
nucleos = multiprocessing.cpu_count()

# GUI Interface
app = tkinter.Tk()
app.title("")
app.geometry('600x600+100+100')
tab_parent = ttk.Notebook(app)
tab1 = tkinter.Frame(tab_parent)
tab2 = tkinter.Frame(tab_parent)
tab2.config(bg="#F7F7F7")
tab1.config(bg="#F7F7F7")
tab_parent.add(tab1, text="tRNAstudio")
tab_parent.pack(expand=1, fill='both')
lbl=tkinter.Label(app, text="")


# Functions definitions
def download_Genome():
    '''
    This function calls download_genome_index.sh to download the Human Genome (hg38) 
    and builds the index for the reference files.
    '''
    lbl.config(text="Downloading and extracting reference genome. This takes a while and freezes the app, don't close it!")
    print ('Downloading Human Reference Genome (hg38), and building the genome indexes... ')
    os.system('bash download_genome_index.sh')
    
    if "genome.4.bt2" and "precursor_tRNA_refgenome.4.bt2" and "mature_fam_tRNA_refgenome.4.bt2" and "mitochondrial_tRNA_refgenome.4.bt2" in os.listdir(app.sourceFolder+"Reference_Genomes/"):
        mb.showinfo("Message", "The genome was downloaded correctly.")
    else:
        mb.showinfo("Message", "The genome was NOT downloaded correctly!!")

 

def retrieve_SRR():
    '''
    This function dowloads sequencing reads in fastq format from NCBI's Sequence Read 
    Archive (SRA).
    '''

    lbl.config(text="Downloading your selected fastq.gz file.")
    SRR = text_Widget.get("1.0",'end-1c')
    
    print ('Downloading sample: '+SRR)
    os.system("fastq-dump " + SRR + " -O " + app.sourceFolder+"/Fastq_downloaded")
    
    if SRR+".fastq" in os.listdir(app.fastqFolder):
        mb.showinfo("Message", "The sample was downloaded correctly")
        print ("The sample was downloaded correctly.")
    else:
        mb.showinfo("Message", "The sample was not dowloaded correctly, check the ID")
        print ("The sample was not dowloaded correctly, check the ID.")
    
    os.chdir(app.scriptsFolder)


def create_text():
    '''
    This function will open a .txt with the sample information (fastq files in Fastq_dowloaded), 
    user will add sample specifications/meatadata.
    '''
    
    sourceFolder = app.sourceFolder
    os.chdir(sourceFolder+"/Fastq_downloaded")
	
    if not any(fname.endswith('.fastq') for fname in os.listdir(os.getcwd())):
    	mb.showinfo("Message", "First download the fastq files or add the fastq files needed for the analysis on the Fastq_downloaded folder.")
    	print ("First download the fastq files or add the fastq files needed for the analysis on the Fastq_downloaded folder.")
    
    else:
        files = os.listdir(os.getcwd())
        if "sample_data_raw.txt" in files:
            if "Linux" in o_sys or "Darwin" in o_sys:
                os.system("rm sample_data_raw.txt")
            files = os.listdir(os.getcwd())
            files = [x for x in files if not x.startswith("._")]
            files = [x for x in files if ".fa" in x]

            
        sample_data = open("sample_data_raw.txt", "w")
        sample_data.write("ID" + "\t" + "Condition" + "\t" + "PE_SE" + "\t" + "Fwd_Rev" + "\t" + "mergeFileID" + "\n")
        
        for i in range(len(files)-1):
            sample_data.write(files[i][:-6] + "\t" + "\n")
        sample_data.write(files[-1][:-6] + "\t")
        sample_data.close()
        
        if "Linux" in o_sys:
            os.system('gedit sample_data_raw.txt')
        if "Darwin" in o_sys:
            os.system('open -a TextEdit sample_data_raw.txt')

        scriptsFolder=app.scriptsFolder       
        os.chdir(scriptsFolder)

def create_text2():
    '''
    This function will open a .txt with the sample information, 
    user will add sample specifications/metadata.
    '''
    sourceFolder=app.sourceFolder
    os.chdir(sourceFolder+"/Fastq_downloaded")
    if "Linux" in o_sys:
        os.system('gedit sample_data.txt')
    if "Darwin" in o_sys:
        os.system('open -a TextEdit sample_data.txt')

    scriptsFolder=app.scriptsFolder       
    os.chdir(scriptsFolder)



def run_alignment():
    '''
    This function performs the alignment pipeline for each sample.
    '''
    os.chdir(app.fastqFolder)
    time_init = time()

    # Check input file and create sample_data.txt
    check_input_file()
    
    sampleData = open("sample_data.txt", "r")   
    sampleData = sampleData.readlines()
    
    if sampleData[-1] == "\n":
        sampleData = sampleData[0:-1]
    sampleData = sampleData[1:]

    for sample in sampleData:
        sample = sample.split("\t")[0]
        sample = sample + ".fastq"

        if 'Linux' in o_sys or 'Darwin' in o_sys:
            os.chdir(app.scriptsFolder)            
            
            # Alignment Whole Genome
            os.system('python3 Aln_WG.py '+app.fastqFolder+ ' ' +sample)

            # Alignment Mature Genome
            os.system('python3 Aln_MG.py '+ sample)

            # Alignment Precursor Genome
            os.system('python3 Aln_PG.py '+ sample)
            
            # Alignemt Mature Genome tolerating more mismatches. 
            os.system('python3 Aln_M1G.py '+ sample)
            
            print ('Alignments done!'+'\n')
            sample = sample.split('.')[0]
            
            print ('Checking mapping results ...')
            os.chdir(app.scriptsFolder)
            check_mapping(sample)
            
            print ('Obtaining counts...')
            os.chdir(app.scriptsFolder)
            counts(sample)
           
            print ('Analyzing mapping quality (MAPQ)...')
            os.chdir(app.scriptsFolder)
            mapping_quality(sample)
            
            print ('Obtaining pileup format ...')
            os.chdir(app.scriptsFolder)
            pileup(sample)
            
            # Prepare samples for R analysis
            r_preparation(sample)
            
            print ('\n'+'Analysis of sample '+sample+' done !'+'\n')

            #Remove temp files
            remove_alignments(sample)

    os.chdir(app.scriptsFolder)
    print("Execution time of the alignments: ", (time() - time_init) / 60, "minutes." )
    mb.showinfo("Message", "Alignment done!")



def check_input_file():
    '''
    This function will check the .txt file containing the metadata,
    and create a new file with an specific format for data processing (unic ids
    for PE data...). It also includes the merging PE fastq file.
    '''
    
    sample_data = open("sample_data_raw.txt", "r")
    sampleData = sample_data.readlines()
    os.chdir(app.fastqFolder)

    if sampleData[-1] == "\n":
        sampleData = sampleData[0:-1]
    sampleData = sampleData[1:]

    NewSampleData= open("sample_data.txt", "w")
    NewSampleData.write("ID" + "\t" + "Condition" + "\n")
    PEsamplesDone = []

    for sample_info in sampleData:
        PE_SE = sample_info.split("\t")[2].replace("\n","")
        condition = sample_info.split("\t")[1].replace("\n","")
        
        # Paired Data
        if PE_SE == "PE":
            sample = sample_info.split("\t")[4].replace("\n","")
            PE_samples = [x for x in sampleData if sample in x]
            if sample not in PEsamplesDone:
                for PE_sample in PE_samples:
                    if 'Fwd' in PE_sample:
                        sampleFwd = PE_sample.split("\t")[0].replace("\t","") + ".fastq"
                    elif 'Rev' in PE_sample:
                        sampleRev = PE_sample.split("\t")[0].replace("\t","") + ".fastq"
                    
                NewSampleData.write(sample.split('.')[0] + "\t" + condition + "\n")
                PEsamplesDone.append(sample)
                
                # Merge PE reads with PEAR
                print ('0-Merging Paired-End reads: '+sample)
                os.system('pear -n 10 -j '+str(nucleos)+' -f '+sampleFwd+' -r '+sampleRev+' -o '+sample)
                os.system('rm -rf *unassembled* *discarded*')
                os.system('mv '+sample+'.assembled.fastq '+sample+'.fastq')
                
        # Single end 
        elif PE_SE == "SE":
            sample = sample_info.split("\t")[0].replace("\t","") + ".fastq"
            NewSampleData.write(sample.split('.')[0] + "\t" + condition + "\n")
        else:
            mb.showinfo("Message", "Pleas fill up the PE_SE column in order to indicate if the sample is Paired-End or Single-End.")




def check_mapping(sample_name):
    
    '''
    Function to summaryze general mapping results
    '''

    sample_data = open("../Fastq_downloaded/sample_data.txt", "r").readlines()
    
    for i in range(len(sample_data)):
        if sample_data[i].split("\t")[0] == sample_name:
            condition = sample_data[i].split("\t")[1].strip("\n")
            break
    
    os.chdir('../Results/'+ sample_name)
        
    readsMapped = int(subprocess.check_output('samtools view -c Alignment_WG/'+sample_name+'_WGloc_mapped_sort.bam', shell=True).decode("utf-8").strip("\n"))
    readsProcessed = int(subprocess.check_output('samtools view -c Alignments/'+sample_name+'_processed.bam', shell=True).decode("utf-8").strip("\n"))
    readsPrecursor = int(subprocess.check_output('samtools view -c Alignments/'+sample_name+'_precursor_sort.bam', shell=True).decode("utf-8").strip("\n"))
    readstRNA = readsPrecursor+readsProcessed
    
    if not os.path.exists('../R_files/'):
        os.makedirs('../R_files/') 

    if not os.path.isfile('../R_files/Alignment_information.txt'):
        fout = open("../R_files/Alignment_information.txt","a")
        fout.write('ID\tCondition\tTotal Reads Mapped\tReads Mapped to tRNA\t% total tRNAs read\t%Processed tRNA reads\n')
    
    fout = open("../R_files/Alignment_information.txt","a")
    
    fout.write(str(sample_name) + "\t" + str(condition) + "\t" + str(readsMapped) +
               "\t" + str(readstRNA) + "\t" + str(format(readstRNA/readsMapped, ".2f")) +
               "\t" +  str(format(readsProcessed/readstRNA, ".2f")) + "\n")



def remove_alignments(sample):
    """
    This function removes the temporary (tmp) aligment files.
    """
    sample = sample.split('.')[0]
    os.chdir("../Results")
    if "Linux" in o_sys or "Darwin" in o_sys:
        os.system("rm -r " + sample+'/Alignment_WG')
        os.system("rm -r " + sample+'/Alignment_MG')
        os.system("rm -r " + sample+'/Alignment_M1G')
        os.system("rm -r " + sample+'/Alignment_PG')
    os.chdir(app.scriptsFolder)




def r_preparation(sample):
    '''
    Files setup for R analysis
    '''

    sample=sample.split('.')[0]
    os.chdir(app.scriptsFolder)
    os.chdir("..")
    if "R_files" not in os.listdir(os.getcwd()+"/Results"):    
        os.mkdir(os.getcwd()+"/Results/R_files")
    
    os.chdir(os.getcwd()+"/Results/R_files")
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
    file = open("../../Fastq_downloaded/sample_data.txt", "r")
    if file.readlines()[-1] != "\n":
        file.close()
        file = open("../../Fastq_downloaded/sample_data.txt", "a")    
        file.write("\n")
        file.close()
    file.close()
    samples = samples[1:]  # Remove the header.

    if "Linux" in o_sys or "Darwin" in o_sys:
        os.chdir("..")
        
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
        os.chdir(app.scriptsFolder)




def r_analysis():
    '''
    This function call the count and modification analysis Rscripts. And launches an Rscript 
    for creating heatmaps depending on the modification ratio and the final report.
    '''
    
    file = open("../Fastq_downloaded/sample_data.txt", "r")
    if file.readlines()[-1] != "\n":
        file.close()
        file = open("../Fastq_downloaded/sample_data.txt", "a")    
        file.write("\n")
        file.close()
    file.close()

    data = pd.read_csv('../Fastq_downloaded/sample_data.txt', sep = '\t')
    print ("Running data analysis: "+'\n')
    print (data)

    subprocess.call(["Rscript", "--vanilla", "Join_Results.r"])
    subprocess.call(["Rscript", "--vanilla", "General_Plots.r"])
    subprocess.call(["Rscript", "--vanilla", "CountsTotal_Plots.r"])
    subprocess.call(["Rscript", "--vanilla", "CountsProcessedPrecursor_Plots.r"])
    subprocess.call(["Rscript", "--vanilla", "ModificationsCoverageAnalysis_Plots.r"])
    subprocess.call(["Rscript", "--vanilla", "ModificationsAnalysis_Heatmaps.r"])
    subprocess.call(["Rscript", "--vanilla", "ModificationsComparison.r"])
    subprocess.call(["Rscript", "--vanilla", "DGE_Analysis.r"])

    mb.showinfo("Message", "Analysis done! You can check the results in the results folder.")



# WIDGETS FOR THE GUI
app.scriptsFolder = os.getcwd()
app.sourceFolder = app.scriptsFolder[:-7]

os.chdir(app.sourceFolder)

try:
    os.mkdir("Fastq_downloaded")
except:
    pass

app.fastqFolder = app.sourceFolder + "Fastq_downloaded"
os.chdir(app.scriptsFolder)

myFont = tkFont.Font(family="Calibri", size=15)
myFont2 = tkFont.Font(family="Calibri", size=14)
myFont3 = tkFont.Font(family="Calibri", size=12)
myFont4 = tkFont.Font(family="Calibri", size=13)
myFont5= tkFont.Font(family="Calibri", size=11)
myFont6= tkFont.Font(family="Calibri", size=10)



if "Darwin" in o_sys:
    text1=tkinter.Label(text="SET UP", width=10,height=1, bg="#F7F7F7")
    text1.place(x=110, y=60)
    text1.config(font='Calibri 14 bold')

    programs_Button=tkinter.Button(tab1, text='Download'+'\n' +'Human Genome (hg38)', width=17, height=2,command=download_Genome, highlightbackground='#a9dbd2')
    programs_Button.place(x=110, y=50)
    programs_Button.config(font=myFont3)
    
    text1=tkinter.Label(text="", width=15, height=1,bg="#F7F7F7")
    text1.place(x=305, y=60)
    text1.config(font='Calibri 14 bold')

    text_Widget = tkinter.Text(height=1, width=13)

    text_Widget.config(font=myFont3)
    text_Widget.insert(tkinter.END, "e.g ERR705691")
    text_Widget.pack()
    def default(event):
        current = text_Widget.get("1.0", tkinter.END)
        if current == "e.g ERR705691\n":
            text_Widget.delete("1.0", tkinter.END)
        elif current == "\n":
            text_Widget.insert("1.0", "e.g ERR705691")

    text_Widget.bind("<FocusIn>", default)
    text_Widget.bind("<FocusOut>", default)
    text_Widget.place(x=330, y=145)

    download_Button=tkinter.Button(tab1, text='Download Sample', width=14, height=2, command=retrieve_SRR, highlightbackground='#a9dbd2')
    download_Button.place(x=300, y=50)
    download_Button.config(font=myFont3)

    text2=tkinter.Label(text='Sample Run Accesion ID', width=16, height=1,bg="#F7F7F7")
    text2.place(x=325, y=130)
    text2.config(font='Calibri 10 italic')

    text3=tkinter.Label(text="tRNAstudio PIPELINE", width=15, height=1,bg="#F7F7F7")
    text3.place(x=240, y=240)
    text3.config(font='Calibri 14 bold')

    text4=tkinter.Label(text='STEP 1', width=5, height=1,bg="#F7F7F7")
    text4.place(x=135, y=260)
    text4.config(font='Calibri 12 bold')

    submit_button = tkinter.Button(tab1, text="Run tRNA Alignments", width = 16, height = 2, command=run_alignment, highlightbackground='#6899aa')
    submit_button.place(x=280, y=245)
    submit_button.config(font=myFont2)
    
    open_text = tkinter.Button(tab1, text='Select and label samples'+'\n'+'for tRNA alignment', width=17, height=2,command=create_text,highlightbackground='#a9dbd2')
    open_text.place(x=110, y = 250)
    open_text.config(font=myFont5)
    
    text1=tkinter.Label(text='STEP 2', width=5, height=1,bg="#F7F7F7")
    text1.place(x=135, y=350)
    text1.config(font='Calibri 12 bold')

    b_countsAnalysis=tkinter.Button(tab1, text='Run Data Analysis', width=16, height=2, command=r_analysis, highlightbackground='#6899aa')
    b_countsAnalysis.place(x=280, y= 335)
    b_countsAnalysis.width=250
    b_countsAnalysis.config(font=myFont2)

    open_text = tkinter.Button(tab1, text='Select samples'+'\n'+'for data analysis', width=17, height=2,command=create_text2,highlightbackground='#a9dbd2')
    open_text.place(x=110, y = 340)
    open_text.config(font=myFont5)

    text1=tkinter.Label(text="Gene Translation Laboratory | IRB Barcelona", width=40, height=1,bg="#F7F7F7")
    text1.place(x=190, y=530)
    text1.config(font='Calibri 9 italic')

    text1=tkinter.Label(text="Integrated pipeline for the analisis of tRNA-Seq Datasets", width=50, height=1,bg="#F7F7F7")
    text1.place(x=160, y=550)
    text1.config(font='Calibri 9 italic')
    app.mainloop()

if "Linux" in o_sys:
	text1=tkinter.Label(text="Set Up", width=10,height=1, bg="#F7F7F7")
	text1.place(x=60, y=45)
	text1.config(font='Calibri 11 bold')

	genome_Button=tkinter.Button(tab1, text='Download'+'\n' +'Human Genome (hg38)', width=17, height=2, command=download_Genome, highlightbackground='#a9dbd2',highlightthickness= 3)
	genome_Button.place(x=80, y=50)
	genome_Button.config(font=myFont5)
	
	genome_Button=tkinter.Button(tab1, text='Download Sample', width=13, height=2, command=retrieve_SRR, highlightbackground='#a9dbd2',highlightthickness= 3)
	genome_Button.place(x=360, y=50)
	genome_Button.config(font=myFont5)


	text_Widget = tkinter.Text(height=1, width=13)
	text_Widget.config(font=myFont5)
	text_Widget.insert(tkinter.END, "e.g ERR705691")
	text_Widget.pack()
	def default(event):
		current = text_Widget.get("1.0", tkinter.END)
		if current == "e.g ERR705691\n":
			text_Widget.delete("1.0", tkinter.END)
		elif current == "\n":
			text_Widget.insert("1.0", "e.g ERR705691")

	text_Widget.bind("<FocusIn>", default)
	text_Widget.bind("<FocusOut>", default)
	text_Widget.place(x=370, y=135)

	text2=tkinter.Label(text='Sample Run Accesion ID', width=23, height=1,bg="#F7F7F7")
	text2.place(x=350, y=120)
	text2.config(font='Calibri 10')


	text3=tkinter.Label(text="tRNAstudio PIPELINE", width=20, height=1,bg="#F7F7F7")
	text3.place(x=230, y=240)
	text3.config(font='Calibri 11 bold')

	text4=tkinter.Label(text='STEP 1', width=6, height=1,bg="#F7F7F7")
	text4.place(x=100, y=270)
	text4.config(font='Calibri 11 bold')

	submit_button = tkinter.Button(tab1, text="Run tRNA Alignments", width = 17, height = 2, command=run_alignment, highlightbackground='#6899aa',highlightthickness= 3)
	submit_button.place(x=310, y=270)
	submit_button.config(font=myFont5)
	
	open_text = tkinter.Button(tab1, text='Select and label samples'+'\n'+'for tRNA alignment', width=22, height=2,command=create_text,highlightbackground='#a9dbd2',highlightthickness= 3)
	open_text.place(x=100, y = 272)
	open_text.config(font=myFont6)
	
	text1=tkinter.Label(text='STEP 2', width=6, height=1,bg="#F7F7F7")
	text1.place(x=100, y=350)
	text1.config(font='Calibri 11 bold')

	b_countsAnalysis=tkinter.Button(tab1, text='Run Data Analysis', width=17, height=2, command=r_analysis, highlightbackground='#6899aa',highlightthickness= 3)
	b_countsAnalysis.place(x=310, y= 350)
	b_countsAnalysis.width=250
	b_countsAnalysis.config(font=myFont5)

	open_text = tkinter.Button(tab1, text='Select samples'+'\n'+'for data analysis', width=22, height=2,command=create_text2,highlightbackground='#a9dbd2',highlightthickness= 3)
	open_text.place(x=100, y = 352)
	open_text.config(font=myFont6)

	text1=tkinter.Label(text="Gene Translation Laboratory | IRB Barcelona", width=45, height=1,bg="#F7F7F7")
	text1.place(x=140, y=530)
	text1.config(font='Calibri 9')

	text1=tkinter.Label(text="Integrated pipeline for the analisis of tRNA-Seq Datasets", width=60, height=1,bg="#F7F7F7")
	text1.place(x=90, y=550)
	text1.config(font='Calibri 9')
	app.mainloop()
