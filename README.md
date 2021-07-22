# tRNAstudio

tRNAstudio is a tool designed to analyze small RNA-seq datasets (Illumina single-end) in order to characterize tRNAs. Includes a specific mapping workflow and provides a report that contains information about tRNA quantification, classification (pre-tRNAs and processed tRNAs), sequence coverage, post-transcriptional tRNA modification levels and differentially expressed genes between two conditions/groups (DESeq2 and Iso-tRNA-cp).

The pipeline is integrated into a freely-available graphical user interface (GUI).

## REQUIREMENTS

**Operative system:** Mac (OS X El Capitan or higher) or Linux.  

**Depencencies**:  
    Python 3.x  
    R 4.1.0  
    Conda  
    bowtie2  
    samtools  
    bedtools  
    gedit/TextEdit  
    picard  
    wget  
    Python packages: pysamstats, pysam, tkinter, pandas , biopython  


## HOW TO USE IT.

**Download the GitHub reposiroty (tRNA-Pipeline):**
1.  Clone the repository or download the zip file. 
2. Unzip the downloaded folder and open the main file. 
3. Open a terminal inside the Scripts folder (right button of the mouse --> open in a terminal). You need to always run the pipeline from the Scripts folder !

**Download dependencies:**

Info: You have to do this step only once. 

You don't have to worry about downloading the required dependencies one by one you can run this command on the terminal in order to install all the dependencies.

` bash requirements.sh`

It first installs Conda that requires downloading the install file from internet and when it is done the user can read the Anaconda User License Agreement (or skip it) pressing ENTER. When it finished the user must type yes in response of "Do you accept the license terms? [yes|no], and press ENTER to confirm the location of the downloaded files. To the question "Do you wish the installer to initialize Anaconda3 by running conda init? [yes|no]" answer no. Please enter the password and type "y" when the terminal ask for it.

Then the script will create an enviroment inside Conda with all the required software and packages (python, R , bowtie2 ....). It is necesary to be aware because it is possible that the software asks for user password and user confirmation in order to install the modules. Then, if there aren't any problems with the python packages you are ready to use the GUI!  


**Use the GUI:**  
Each time that you want to use the GUI you need to first activate the conda enviroment (it contains all the software and programs needed in order to run the tRNA pipeline). Copy the next command on the terminal:

`conda activate tRNAstudioEnv` or `source activate tRNAstudioEnv`

Then run the GUI:

`python GUI.py`

Once the GUI is working you will have the next options:

**Download Human Genome.** Downloads human genome Hg38 from USCS, unzip and index it, all the reference files will be saved in the "Reference_Genomes" folder. This process takes some time, but when it's done the user will be notified. Not user interaction needed.

**Download sample.** This button download the .fastq file from the European Nucleotide Archive (ENA Browser) from its Run Accession ID. Run Accession ID examples: SRR7216347 or ERR705691. First it is necessary to enter the sample ID (without blanks or tabulators) and then press the button. Once the .fastq file is downloaded and unziped in the "Fastq_downloaded" a message will notify the user, and then more samples can be downloaded following the same process. In the case that you have the .fastq files in your computer you don't need to download them from the EBI, just copy paste them in the "Fastq_downloaded" folder.

**Select and label samples for tRNA alignment**  Indicate sample information. This button opens a .txt file with a text editor that contains the sample's ID and a column (Condition) in blank that MUST be filled by the user indicating the group (for example: Control/Treated) to which each sample belongs. In each row the information of the 2 columns must be separated by 1 tabulator (\t) and without blank spaces. Don't introduce any blank in the file, and finally save changes and close the file.

**Run tRNA Alignments** Pressing this button it is perfomed the several rounds of alignment for each sample against the whole human genome, the mature tRNA genome and the precursor tRNA genome. This process needs a lot of computational power so the user is asked not to do another important thing with the computer while this process is running. This process can last several hours and it requires a lot of computational resources so it is recommended to the user not to use the computer while this step is runing, and finally will notify the user when it is finished.

**Select samples for data analysis** This button opens the previously saved .txt file with the samples information, allowing the user to edit the file and select which samples to analyze. 
  
**Run Data Analysis** Once all the alignments are done the used has to press this button to compute the counts analysis, modification analysis for each group and the differential gene expression analysis between the groups. Finally it generates a report per group and a report with the comparison between groups summarizing the most important information. This step takes time too, and the user it's informed when it finished.

<p align="center">
  <img src="https://github.com/GeneTranslationLab-IRB/tRNAstudio/blob/main/tRNAstudioGUI.jpg" />
</p>

## RESULTS

In the results folder are generated the following folders:
Counts Plots.

In this folder are contained the count plots:
- Total counts by isodecoder tRNA and by tRNA gene grouped by isoacceptor tRNA in each condition.
- Proportion mature/precursor by isodecoder tRNA and by tRNA gene group by isoacceptor tRNA in each condition.
DEG.

Differential gene expression. In this folder there are results of the differential gene expression analysis like the heatmaps with the gene expression between groups and files with DESeq2 and iso-tRNA-CP results.
Heatmaps.

Heatmaps produced by "heatmaply" [1] with the information of the modification ratio in each position of every tRNA gene.
Modification ratio plots.

Plots showing the coverage and the modification ratio by position in each tRNA gene. In "By_family" folder there are a summary of the coverage by tRNA isoacceptor.
Plots.

Folder with general plots as the number of reads by sample, PCA, a barplot with the reads with mapping quality higher and lower than 2 MAPQ. isotRNACP folder contains plots with proportion of each gene against its isodecoder pool.
Reports.

In this folder are contained the reports of the pipeline, a report with the characterization for each group and one report with the comparison of the groups.
R_files.

Internal files used by Rscripts to compute our analysis.
Folders with the ID of each sample.

It contains the results of the alignment of each sample, the number of counts and the base calling.
Support

If you have any questions or issues, please use Issues Section
References.

    Galili, Tal, O'Callaghan, Alan, Sidi, Jonathan, Sievert, Carson (2017). “heatmaply: an R package for creating interactive cluster heatmaps for online publishing.” Bioinformatics.
