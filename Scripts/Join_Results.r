#! /usr/bin/Rscript

suppressPackageStartupMessages({
library(dplyr)
})

sample_data = read.delim("../Fastq_downloaded/sample_data.txt", header=T)

setwd("../Results/")
first_name = paste0(sample_data[1,1],"/Counts/",
                    sample_data[1,1],"_counts_total.txt" )
first_file = read.table(first_name, header = T)
first_file$counts = apply(first_file[,c(2,3)], 1, function(x){ sum(x)}  )
first_file = data.frame(first_file$tRNA_ID, first_file$counts)
data_merged = first_file
names_files = c(sample_data[1,1])

for(i in c(2:nrow(sample_data))){
  name_file = paste0(sample_data[i,1], "/Counts/",sample_data[i,1],
                     "_counts_total.txt")
  names_files = c(names_files, sample_data[i,1])
  data_file = read.table(name_file, header = T)
  data_file$counts = apply(data_file[,c(2,3)], 1, function(x){ sum(x)}  )
  data_merged = cbind(data_merged, as.numeric(data_file$counts))
}

names_genes = data_merged[,1]
data_merged = data.frame(data_merged[,2:ncol(data_merged)])
colnames(data_merged) = c(names_files)
rownames(data_merged) = names_genes
dir = "R_files/Counts/"
write.table(data_merged, file=paste0(dir,"Counts_by_gene_total_for_DESeq2.txt"), 
            row.names = T, sep="\t", quote=F)


setwd("../Scripts")
