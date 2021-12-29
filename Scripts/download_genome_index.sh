#!/bin/bash

#Changing the directory to download to Reference_Genomes

cd ..
mkdir Reference_Genomes
cd Reference_Genomes
#Downloading the file from illumina and uncompressing it
wget http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Homo_sapiens/UCSC/hg38/Homo_sapiens_UCSC_hg38.tar.gz
tar -zxvf Homo_sapiens_UCSC_hg38.tar.gz
#Moving all the needed files from their folder to the main Reference_Genomes folder
mv Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome.1.bt2 genome.1.bt2
mv Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome.2.bt2 genome.2.bt2
mv Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome.3.bt2 genome.3.bt2
mv Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome.4.bt2 genome.4.bt2
mv Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa genome.fa
mv Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.dict genome.dict
mv Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa.fai genome.fa.fai
mv Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/GenomeSize.xml GenomeSize.xml
mv Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome.fa genome.fa
mv Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome.rev.1.bt2 genome.rev.1.bt2
mv Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome.rev.2.bt2 genome.rev.2.bt2
#Deleting the tar.gz file
rm Homo_sapiens_UCSC_hg38.tar.gz
#Deleting the rest
rm -r Homo_sapiens
#Building indexs of the "genome" from tRNA families and precursor tRNA families
bowtie2-build mature_fam_tRNA_refgenome.fa mature_fam_tRNA_refgenome
bowtie2-build precursor_fam_tRNA_refgenome.fa precursor_fam_tRNA_refgenome
bowtie2-build mitochondrial_tRNA_refgenome.fa mitochondrial_tRNA_refgenome