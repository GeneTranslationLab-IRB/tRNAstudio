# tRNAstudio

tRNAstudio is a tool designed to analyze small RNA-seq datasets (Illumina single-end) in order to characterize tRNAs. Includes a specific mapping workflow and provides a report that contains information about tRNA quantification, classification (pre-tRNAs and processed tRNAs), sequence coverage, post-transcriptional tRNA modification levels and differentially expressed genes between two conditions/groups (DESeq2 and Iso-tRNA-cp).

The pipeline is integrated into a freely-available graphical user interface (GUI)

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

## EXPLANATION OF THE GUI

**Download Genome.** Downloads human genome Hg38 from USCS, unzip and index it. This process takes some time, but when it's done the user will be notified.

**Download sample.** This button download the .fastq file from the European Bioinformatic Institute (EBI) from its Run Accession ID. Run Accession ID examples: SRR7216347 or ERR705691. First it is necessary enter the ID and then press the button. Once the sample is downloaded and unzipped a popped message will notify it and it is possible to download more samples. In the case that you have the .fastq files in your computer you don't need to download them from the EBI, just copy paste them in the "Fastq_downloaded" folder.

**Indicate sample information.**  Indicate sample information. This button opens a .txt file with a text editor that contains the sample's ID and a column (Condition) in blank that MUST be filled by the user indicating the group in which each sample belongs. In each row the information of the 2 columns must be separated by 1 tabulator (\t) and without blank spaces.

  **Run aligment.** Pressing this button it is perfomed the several rounds of alignment for each sample against the whole human genome, the mature tRNA genome and the precursor tRNA genome. This process needs a lot of computational power so the user is asked not to do another important thing with the computer while this process is running. This process can last several hours and will notify the user when it is finished.
  
**Results Report.** Perform the counts analysis, modification analysis for each group and the differential gene expression analysis between the groups. Finally it generates a report per group and a report with the comparison between groups summarizing the most important information.

## HOW TO USE IT.

**Download the GitHub reposiroty (tRNA-Pipeline):**
1.  Clone the repository or download the zip file. 
2. Unzip the downloaded folder and open the main file. 
3. Open a terminal inside the Scripts folder (right button of the mouse --> open in a terminal). You need to always run the pipeline from the Scripts folder !

**Download dependencies:**
Info: You have to do this step only once. 
You don't have to worry about downloading the required dependencies one by one you can run this command on the terminal in order to install all the dependencies.

` bash requirements.sh`

It first installs conda (open source package management system and environment management system) automatically, that requires downloading the install file from internet and when it is done the user can read the Anaconda User License Agreement (or skip it) pressing ENTER. When it finished the user must type yes in response of "Do you accept the license terms? [yes|no], and press ENTER to confirm the location of the downloaded files. To the question "Do you wish the installer to initialize Anaconda3 by running conda init? [yes|no]" answer no. Please enter the password and type "y" when the terminal ask for it.

Then the script will create an enviroment inside conda with all the required software and packages (python, R , bowtie2 ....). It is necesary to be aware because it is possible that the software asks for user password and user confirmation in order to install the modules. Then, if there aren't any problems with the python packages you are ready to use the GUI!  


**Use the GUI**  
Each time that you want to use the GUI you need to first activate the conda enviroment (it contains all the software and programs needed in order to run the tRNA pipeline). Copy the next command on the terminal:

`conda activate tRNA-PipelineEnv`

Then run the GUI:

`python GUI.py`

Once the GUI is working you will have the next options:

**Download Genome:**
With this button the Human Genome 38 will be downloaded, unziped and indexed in the Reference_Genomes folder. Not user interaction needed.

**Download sample:**
Type the Run Accession ID from the EBI (SRRXXXXXXX) in the box (without blanks or tabulators) and press the button Download sample. The fastq file will be downloaded and unziped in the Fastq_downloaded folder. This process must be repeated until all the samples of the experiment had been downloaded.

**Indicate sample information:**
Pressing this button it is opened the file sample_data.txt located in the "Fastq_downloaded" folder. This file consists in a table of 2 columns, the first one with the ID of each sample and in the second column the user MUST write the group in which each sample belongs (for example: Control/Treated). Don't introduce any blank in the file, and finally save changes and close the file.

**Run alignment**
This button performs the several rounds of alignent of each sample agains the whole human genome, the precursor tRNA genome and the mature tRNA genome. The user don't need to press anything but the button. This step can last several hours and it requires a lot of computational resources so it is asked to the user not to use the computer while this step is runing.

**Results Report:** 
Once all the alignments are done press this button to compute the analysis of counts, modification ratio and differential gene expression. This step takes time too, and the user it's informed when it finished.

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
