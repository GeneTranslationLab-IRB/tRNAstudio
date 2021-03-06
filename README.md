<div align="center">
    <h1 align="center">tRNAstudio</h3>
  </p>
</div>

tRNAstudio is a tool designed to analyze small RNA-seq datasets (single-end or paired-end data) in order to characterize tRNAs. It includes a specific mapping workflow and extracts information on tRNA gene expression (DESeq2 and Iso-tRNA-CP), classification of sequencing reads likely derived from precursor tRNAs (pre-tRNAs) or mature tRNAs (processed tRNAs), tRNA gene sequence coverage, and post-transcriptional tRNA modification levels.
<br /> The pipeline is integrated into a freely-available graphical user interface (GUI) that simplifies the analyses for users without programming skills.

## REQUIREMENTS

**Operating system:** Mac (OS X El Capitan or higher) or Linux.  



## HOW TO USE IT

**Download the GitHub repository (tRNAstudio):**
1. Click on the green button "code" at the top of the page and download the zip file.
2. Unzip the downloaded folder and open the "tRNAstudio-main" folder. (The folder can be moved to another directory, make sure that the new directory/location does not contain a folder with a name that has blank spaces)
3. For Linux Users: Click the "Scripts" folder with the right button and select "Open in a terminal". 
<br /> For Mac Users: Click the "Scripts" folder with the right button followed by “Services” and select “New Terminal at Folder.”
<br /> <br /> <sup> **Note:**  It is always necessary to run the pipeline from the "Scripts" folder!  </sup>

**Download dependencies:** 

This step is only done on the first time the pipeline is being installed. 

Run this command on the terminal.

`bash requirements.sh`

This command automatically installs Conda (an open source package management system). Follow the next steps in order to complete the installation. 

0. This step is only for Mac users: First, it can ask for your password "Checking for sudo access (which may request your password)... Password:". Writhe your password and press "ENTER". Press "ENTER" again to continue.
1. Read and accept the Anaconda User License Agreement on the terminal. Press "ENTER" to go through the license terms until reaching the question "Do you accept the license terms? [yes|no], answer **"yes"** and press "ENTER". 
2. Press "ENTER" to confirm the location of the downloaded files. 
3. To the question "Do you wish the installer to initialize Anaconda3 by running conda init? [yes|no] answer **"no"** and press "ENTER". If needed, Enter your user password and type "y" when the terminal asks for it.
4. The script will then create an environment inside Conda with all the software and packages required to run tRNAstudio. Some of these software may ask for the user password or confirmation in order to install the modules; therefore the user must be attentive to these instances and supply this information if required. 
5. Once the installation is done, close the terminal and reopen it again from the scripts folder as described above. Now you are  ready to use the tRNAstudio GUI!  

**Use the GUI:**  

The next **two commands** are always needed in order to run tRNAstudio (remember to run the commands from the "Scripts" folder).

 **1.Command to activate the environment**

First, activate the conda environment by running the following command on the terminal:

```conda activate tRNAstudioEnv```

If this error appears "CommandNotFoundError: Your shell has not been properly configured to use ‘conda activate’...". Activate the environment using:

```source activate tRNAstudioEnv```

**2.Command to run the GUI** 

Then run the following command to run the GUI:

```python3 tRNAstudioGUI.py```

Once the GUI is running the user will have the following options (Fig.1):

<p align="center">
  <img src="https://github.com/GeneTranslationLab-IRB/tRNAstudio/blob/main/img/tRNAstudioGUI.jpg" width="500" height="530"/>
    <p align="center">Figure 1. tRNAstudio graphical user interface (GUI)</p>
</p>

<sup> **Note:** The pipeline is designed to run one process at a time, once that you click one of the buttons you have to wait until it finishes in order to start a new process.  </sup>

<sup> **Note:** The user can view the progress of each step in the terminal.  </sup>

- **Download Human Genome (hg38).**  This button downloads the human genome Hg38 from UCSC, and builds the indexes saved in the "Reference_Genomes" folder. This process takes some time (can last 2 hours) and a pop-up dialog will inform the user when the process is finished. No user interaction needed. This process will have to be done only once (provided that the user does not manually delete this downloaded file).

- **Download sample.** This button downloads the .fastq files from the Gene Expression Omnibus (GEO) using its Run Accession ID. Run Accession IDs are those with the prefixes SRR…, ERR…, or DRR… (e.g.ERR705691). Enter the Run Accession ID (without blank spaces or tabulators) and then press the "Download sample" button. Once the .fastq file is downloaded, it will be automatically unzipped and stored in the "Fastq_downloaded'' folder. A message will notify the user that this process was correctly achieved. The user can then download additional samples (i.e. other Run Accession IDs) following the same procedure. <br /> <sup> **Note:** If the user already has the desired unziped .fastq files in their computer (e.g. their own in-house datasets or datasets obtained upon manual downloading from other repositories), these files should be placed inside the "Fastq_downloaded" folder and will become accessible to be analyzed by tRNAstudio. </sup>

- **Select and label samples for tRNA alignment.**  This button opens a .txt file with a text editor that shows the sample's ID contained in the "Fastq_downloaded'' folder (**Column 1**) and empty columns that MUST be filled by the user (columns might be misaligned given that in the .txt file they are defined by tabulators): <br />
<br /> &nbsp;&nbsp;&nbsp;**Column 2**: Condition. Indicate the group to which each sample belongs (e.g. Control or Treated). Do not use only numerical values to indicate the condition.
<br /> &nbsp;&nbsp;&nbsp;**Column 3**: PE_SE. Specify if the sample is paired-end (write down "PE") or single-end (write down "SE").
<br /> &nbsp;&nbsp;&nbsp;**Column 4**: Fwd_Rev. If the sample is PE, the user will have two files per each sample downloaded. One of these samples will 
<br /> correspond to the Forward R1 and the other one to the Reverse R2. The user has to indicate in this column which file is the Forward R1 (write down "Fwd") and which one is the Reverse R2 (write down "Rev"). 
<br /> &nbsp;&nbsp;&nbsp;**Column 5**: mergeFileID. If the sample is PE, Forward and Reverse files will be automatically merged into one file before running the pipeline. Therefore, in this case the user will have to select a unique name for the merged files (without blank spaces). In this column, the same chosen name for the merged file has to be indicated in both Forward and Reverse files. <br /> <sup> **Note:** if the dataset is SE, the user does not have to fill columns 4 and 5 of the .txt file. </sup> <br /> <sup> **Note:** In each row of the .txt file (Fig. 2), the information of the columns must not contain blank spaces, and must be separated by 1 tabulator (i.e. write down the information and click "tab" to move to the next column in the row). </sup>
<p align="center">
<br /> <img src="https://github.com/GeneTranslationLab-IRB/tRNAstudio/blob/main/img/sample_data.jpg" alt="" width="360" height="260"/> </p>
    <p align="center">Figure 2. File containing sample data</p>
</p>
<br /> <br /> Once finishing providing the required information in each column save the changes (Control + S) and close the file. DO NOT press "ENTER" after completing the last column of the file as this will incorporate an additional blank row in the .txt file that will be interpreted as an error. 


- **Run tRNA Alignments.** This button executes the alignment pipeline implemented in tRNAstudio. This process can last several hours and requires a lot of computational power so we recommend not performing other demanding processes while the pipeline is running. The user will be notified with a pop-up dialogue when this process is finished.

- **Select samples for data analysis.** This button opens the previously saved .txt file with the sample information, allowing the user to select which samples to analyze (e.g. delete samples that are not required for  data analysis) and to indicate the group to which each sample belongs (e.g. Control or Treated). Do not use only numerical values to indicate the condition. If the file is modified, save changes (Control + S) and close it.
 
   <sup> **Note:** Remember that tRNAstudio is designed to analyze and compare two conditions. </sup>

  
- **Run Data Analysis.** This button computes all the parameters that tRNAstudio can assess (e.g. tRNA quantification, modification analysis for each group, differential gene expression analysis between the different conditions, etc.). Press this button after the alignments are done and after the samples that are to be analyzed have been selected. This step is also time-consuming. The user will be notified with a pop-up dialogue when this process is finished.


## RESULTS
<sup> **Note:** We recommend to save the next folders: **General_Plots, Counts_Plots, Modification_Coverage, DGE and Reports**, on another directory named with a specific project name. Do not move **R_files or Folders named after the ID of each sample** from the Results folder </sup>

The following folders are generated inside the "Results" folder:

- **Folders named after the ID of each sample.**

   We recommend to browse this folder in order to perform custom analysis.
    
    This folder contains 3 subfolder:

    - "Alignments": .bam files that contain sequence alignment data for processed, precursor and mitochondrial tRNAs. It also contains sorted and indexed .bam files needed in the case that the user wants to use visualization tools such as the Integrative Genomics Viewer (IGV).

    - "Counts": .txt files with the counts for processed, precursor and mitochondrial tRNAs. Each .txt file contains 3 columns [tRNA ID] [counts with a MAPQ > 2] [counts with a MAPQ ≤ 2].
    - "Base_Calls": .txt file with the number of reads with a given base at each tRNA position (base calling) for each tRNA gene.

- **R_files.**
Internal files used by tRNAstudio to compute the analysis. This folder should be managed only by experienced bioinformaticians. 

- **General_Plots.**
Folder with barplots for the total number of reads mapped to tRNAs in all conditions and separated by condition, proportion of the total number of reads mapped to tRNAs by mapping quality in all conditions and separated by condition and Principal Component Analysis (PCA) plot. 

- **Counts_Plots.**
    This folder contains 3 subfolders:
    - "Total": plots showing the total number of reads mapped to tRNAs by isoacceptor, isodecoder and tRNA gene family for each condition.
    - "Processed_vs_precursor": plots showing the proportion of reads derived from processed tRNA and pre-tRNA by isoacceptor, isodecoder and tRNA gene family  for each condition.
    - "Proportion_by_isodecoder": plots showing the proportion of reads mapping to each tRNA gene family within their isodecoder set.

- **Modification_Coverage.**
    This folder contains .txt files and plots for the analysis of tRNA modifications and sequence coverage. It contains 4 subfolders:
    - "Modification_Coverage_Plots": line charts showing the sequence coverage and the modification ratio by position for each tRNA gene family. In addition, the .txt files with the data used to generate the plots is included in this folder. 
    - "Modification_Heatmaps_Plots": heatmaps showing the modification ratio at selected positions for each tRNA gene family.
    - "Modification_Comparison_Plots": line charts showing the modification ratio by position in each tRNA gene family by condition. Consider this data in the context of gene coverage; as low gene coverage may give misleading results. 
    - "Coverage_Plots": plots showing the coverage for each tRNA gene family (images grouped by isoacceptor).

- **DGE.**
Differential gene expression analysis. 
This folder contains the numerical results of the differential gene expression analysis (e.g. counts used for generating gene expression heatmaps, files with DESeq2 and iso-tRNA-CP results, etc.). It includes images in two different formats; .jpeg and .pdf, for better quality images. An interactive .html summary is provided with the results of DESeq2 and ISO-tRNA-CP.


- **Reports.**
This folder contains a report (interactive .html) with the characterization for each group and one report with the comparison between the two groups. Inside the report you have links (Titles in blue) with additional information on some of the analysis.  

## ADAPTING tRNAstudio TO ANALYZE DATASETS FROM NON-HUMAN SPECIES. 

<sup> **Note:** this section is intended for experienced bioinformaticians. </sup>

tRNAstudio is designed for the analyses of human datasets, however, it can be adapted to analyze datasets from other species. 

First, a database with the sequences of the tRNAs (nuclear and mitochondrial tRNA genes if needed) for the species of interest need to be generated. To do so, tRNA sequences can be downloaded from specific tRNA databases (eg. GtRNAdb, mitotRNAdb). Make sure that the obtained sequences contain the following information:

- tRNA gene name (tRNA ID).
- Sequence.
- Genomic coordinates of the full gene. 
- Coordinates of intronic sequences.

In addition, it will be necessary to indicate the consensus tRNA base position (positions 1 to 73) for each nucleotide on each tRNA gene sequence. This information can be extracted from the structural data associated to each sequence as defined by tRNAscan-SE, and can be validated by browsing dedicated databases (e.g. tRNAviz). 

Next, adapt the files within the “Reference_Genomes” folder with the tRNA data for the species of interest. 

Genomes.
-	**Reference genome:** Download the complete genome for the species of interest from the UCSC Genome Browser.
-	**Precursor genome:**  Generate a .fastq file with Precursor tRNA sequences for each nuclear tRNA gene. To obtain precursor tRNA sequences, use the genomic tRNA gene sequence including 5’-end and 3’-end flanking regions (50 nucleotides upstream and 50 nucleotides downstream of the gene) to simulate leader and trailer regions. Then, group all identical sequences into families so that the generated genome contains unique sequences. To name each tRNA gene/tRNA family (tRNA ID) we recommend using the format “tRNA-AminoAcid-Anticodon-FamilyNumber-GeneNumber” for labeling single genes and “tRNA-AminoAcid-Anticodon-FamilyNumber” for labeling gene families.
-	**Mature genome:** Generate a .fastq file with mature tRNA sequences for each nuclear tRNA gene. Use the tRNA gene sequence without intronic sequences, and with the addition of the trinucleotides “CCA” at their 3’-end (all tRNA genes) and of a “G” at the 5’-end (only for His tRNAs). Identical sequences need to be grouped into families so that the generated genome contains unique sequences. 
-	**Mitochondrial genome:** Generate a .fastq file with each mitochondrial tRNA sequence.

Once all reference genome files have been generated, they need to be indexed using bowtie2-build.

Next, extract the following information (see below) from the generated genomes and create independent reference files with the following names (to create this files use as reference the files contained in the "Reference_Genomes/info" folder):


tRNAsCoordinates.bed: .bed file containing the coordinates for nuclear tRNA gene plus 50 nucleotides flanking regions upstream and downstream of the gene, and mitochondrial tRNA genes.
Columns description: chr| start| end| tRNA ID| score| strand


leader_trailer_intron.bed: .bed file containing the coordinates of leader, trailer and intronic regions for each nuclear tRNA gene. 
Columns description: chr| start| end| tRNA ID |score| strand


intron.bed: .bed file containing the coordinates of intronic regions for each nuclear tRNA gene. 
Columns description: chr| start| end| tRNA ID| score| strand


tRNAsAnnotation.gtf: .gtf file containing the tRNA gene annotation for precursor, mature and mitochondrial tRNA sequences.
Columns description:  |seqname| source| feature |start |end |score |strand| attribute|

tRNAsPositions.txt: .txt file that contains the tRNA ID of each tRNA gene present in the Mature genome, and the corresponding consensus tRNA positions. 

tRNAsID.txt: Dictionary with the precursor ID as key and the mature family ID as value. 

Upon generating all required files, replace the original files from tRNAstudio located in the “Reference_Genomes” and  “Reference_Genomes/info” folders with these new files. tRNAstudio should then be ready to study non-human species; however, we strongly recommend validating the pipeline thoroughly before implementing a modified version of tRNAstudio in an automated and unsupervised manner.

## SUPPORT 

If you have any questions or issues, please use the Issues Section

