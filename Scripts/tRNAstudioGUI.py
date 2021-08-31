


'''
GUI Interface Integrated pipeline for the analisis of tRNA-Seq Datasets (tRNAstudio)
'''


###############################################################################
print ('\n'+'############# tRNAstudio ###############'+'\n')

# Built-in/Generic Imports #

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning) 

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

# Import functions #
from Functions import counts
from Functions import pileup

# System
o_sys=platform.system()


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
    This function calls the scripts to download the Human genome and move 
    the files where they need to be.
    '''
    lbl.config(text="Downloading and extracting reference genome. This takes a while and freezes the app, don't close it!")

    if 'Linux' in o_sys or 'Darwin' in o_sys:
        os.system('bash download_genome_index.sh')
        
    if "genome.4.bt2" and "precursor_tRNA_refgenome.4.bt2" and "families_tRNA_refgenome.4.bt2" in os.listdir(app.sourceFolder+"Reference_Genomes/"):
        mb.showinfo("Message", "The genome was downloaded correctly.")
    else:
        mb.showinfo("Message", "The genome was NOT downloaded correctly!!")

 

def retrieve_SRR():
    '''
    This function recieves an SRR accession number and with string methods 
    creates the ftp url where the file is, and dowloads it. It uses different 
    methods depending on the operating system.
    '''

    lbl.config(text="Downloading your selected fastq.gz file.")
    SRR = text_Widget.get("1.0",'end-1c')
    first_digits=SRR[0:6]
    last_digits='00'+SRR[-1:]
    ftp_link1='ftp://ftp.sra.ebi.ac.uk/vol1/fastq/'+first_digits+'/'+last_digits+'/'+SRR+'/'+SRR+'.fastq.gz'
    ftp_link2 ="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/"+first_digits+"/"+SRR+"/"+SRR+".fastq.gz"

    if 'Linux' in o_sys or 'Darwin' in o_sys:
        os.system('wget -P ' + app.fastqFolder + " " + ftp_link1)
        if SRR+".fastq.gz" not in os.listdir(app.fastqFolder):
            os.system('wget -P ' + app.fastqFolder + " " + ftp_link2)
        
        print ("Unzipping fastq file....")
        os.system("gunzip " + app.fastqFolder + "/" + SRR + ".fastq.gz")
    if SRR+".fastq" in os.listdir(app.fastqFolder):
        mb.showinfo("Message", "The sample was downloaded correctly")
        print ("The sample was downloaded correctly.")
    else:
        mb.showinfo("Message", "The sample was not dowloaded correctly, check the ID")
        print ("The sample was not dowloaded correctly, check the ID.")
    
    os.chdir(app.scriptsFolder)


def create_text():
    '''
    This function will open a txt with the sample information, 
    user will add sample specifications
    '''
    sourceFolder=app.sourceFolder
    os.chdir(sourceFolder+"/Fastq_downloaded")
    files = os.listdir(os.getcwd())
    if "sample_data.txt" in files:
        if "Linux" in o_sys or "Darwin" in o_sys:
            os.system("rm sample_data.txt")
        files = os.listdir(os.getcwd())
        files = [x for x in files if not x.startswith("._")]

        
    sample_data = open("sample_data.txt", "w")
    sample_data.write("ID" + "\t" + "Condition" + "\n")
    for i in range(len(files)-1):
        sample_data.write(files[i][:-6] + "\t" + "\n")
    sample_data.write(files[-1][:-6] + "\t")
    sample_data.close()
    if "Linux" in o_sys:
        os.system('gedit sample_data.txt')
    if "Darwin" in o_sys:
        os.system('open -a TextEdit sample_data.txt')

    scriptsFolder=app.scriptsFolder       
    os.chdir(scriptsFolder)

def create_text2():
    '''
    This function will open a txt with the sample information, 
    user will add sample specifications
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
    This function reads the file with the samples in order to perform the tRNA 
    alignment pipeline.
    '''
    time_init = time()
    os.chdir(app.fastqFolder)
    file = open("sample_data.txt", "r")
    
    sample_data = file.readlines()
   
    sample_data = sample_data[1:]
    for sample in sample_data:
        sample = sample.split("\t")[0]
        sample = sample + ".fastq"
        

        lbl.config(text="Starting the pipeline.")
        lbl.config(text="Done.")
        lbl.config(text="Starting the whole genome alignment. This takes time.")
        if 'Linux' in o_sys or 'Darwin' in o_sys:

            os.chdir(app.scriptsFolder)
            os.system('python3 Aln_WG.py '+sample+ ' ' + app.fastqFolder)
            lbl.config(text="Done!")
            
            lbl.config(text="Aligning versus the mature genome.")
            os.system('python3 Aln_MG.py '+ sample)
            lbl.config(text="Done!")

            lbl.config(text="Aligning versus the precursor genome.")
            os.system('python3 Aln_PG.py '+ sample)
            lbl.config(text="Done!")
            
            lbl.config(text="Aligning versus the mature genome with one mismatch in the seed.")
            os.system('python3 Aln_M1G.py '+ sample)
            lbl.config(text="Done!")
            print ('Alignments done !!')

            os.chdir(app.scriptsFolder)
            lbl.config(text="Obtaining the final counts.")
            sample = sample.split('.')[0]
            counts(sample)
            lbl.config(text="Done!")
            
            os.chdir(app.scriptsFolder)
            lbl.config(text="Doing the pileup.")
            pileup(sample)
            lbl.config(text="Done!")
            r_preparation(sample)

            print ('\n'+'Analysis of sample '+sample+' done !'+'\n')
            
            #remove_alignments(sample)
            
    os.chdir(app.scriptsFolder)
    print("Execution time of the alignments: ", (time() - time_init) / 60, "minutes." )
    mb.showinfo("Message", "Alignment done!")




def remove_alignments(sample):
    """
    This function remove the aligments files.
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
    This script takes all the counts step result files and joins them
    in a single file in order to analyse them in R
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
            "all_total.txt" 
        path_to_file2 = os.getcwd()+"/"+sample+"/Counts/"+sample+\
            "all_mature_sort.txt"                 
        path_to_file3 = os.getcwd()+"/"+sample+"/Base_calling/"+sample+ \
            "_all_base_calling_by_pos_CORRECT_OK.txt" 
        path_to_file4 = os.getcwd()+"/"+sample+"/Alignment_WG/"+sample+\
            "_WGloc_only_trna.bam"
        path_to_file5 = os.getcwd()+"/"+sample+"/Alignment_MG/"+sample+\
            "_MGloc_mapped.bam"
        path_to_file6 = os.getcwd()+"/"+sample+"/Alignment_PG/"+sample+\
            "_PGloc_mapped.bam"
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
    print ("Performing data analysis ")
    
    subprocess.call(["Rscript", "--vanilla", "Join_results.r"])
    
    #Counts analysis.
    subprocess.call(["Rscript", "--vanilla", "PLOTS_Counts_Total.r"])
    subprocess.call(["Rscript", "--vanilla", "PLOTS_Prop_Mature_Precursor.r"])
    
    #Modification analysis.
    subprocess.call(["Rscript", "--vanilla", "PLOTS_Modifications_Analisis.r"])
    
    
    subprocess.call(["Rscript", "--vanilla", "General_plots.r"])
    subprocess.call(["Rscript", "--vanilla", "Heatmaps.r"])
    subprocess.call(["Rscript", "--vanilla", "DEG_analysis.r"])
    
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

    programs_Button=tkinter.Button(tab1, text='Download'+'\n' +'Human Genome (hg38)', width=14 , height=2,command=download_Genome, highlightbackground='#a9dbd2')
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

    download_Button=tkinter.Button(tab1, text='Download Sample', width=10, height=2, command=retrieve_SRR, highlightbackground='#a9dbd2')
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

    submit_button = tkinter.Button(tab1, text="Run tRNA Alignments", width = 14, height = 2, command=run_alignment, highlightbackground='#6899aa')
    submit_button.place(x=280, y=245)
    submit_button.config(font=myFont2)
    
    open_text = tkinter.Button(tab1, text='Select and label samples'+'\n'+'for tRNA alignment', width=14, height=2,command=create_text,highlightbackground='#a9dbd2')
    open_text.place(x=110, y = 250)
    open_text.config(font=myFont5)
    
    text1=tkinter.Label(text='STEP 2', width=5, height=1,bg="#F7F7F7")
    text1.place(x=135, y=350)
    text1.config(font='Calibri 12 bold')

    b_countsAnalysis=tkinter.Button(tab1, text='Run Data Analysis', width=14, height=2, command=r_analysis, highlightbackground='#6899aa')
    b_countsAnalysis.place(x=280, y= 335)
    b_countsAnalysis.width=250
    b_countsAnalysis.config(font=myFont2)

    open_text = tkinter.Button(tab1, text='Select samples'+'\n'+'for data analysis', width=14, height=2,command=create_text2,highlightbackground='#a9dbd2')
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
