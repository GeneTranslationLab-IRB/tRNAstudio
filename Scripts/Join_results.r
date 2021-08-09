
suppressPackageStartupMessages({
library(dplyr)
})

sample_data = read.table("../Fastq_downloaded/sample_data.txt", header=T,
                         stringsAsFactors = F)

setwd("../Results/R_files/Counts")


first_name = paste0(sample_data[1,1],"all_total.txt" )
first_file = read.table(first_name)
data_merged = first_file
names_files = c(sample_data[1,1])

for(i in c(2:nrow(sample_data))){
  name_file = paste0(sample_data[i,1], "all_total.txt")
  names_files = c(names_files, sample_data[i,1])
  data_file = read.table(name_file)
  data_merged = cbind(data_merged, as.numeric(data_file$V2))
}

names_genes = data_merged[,1]
data_merged = data.frame(data_merged[,2:ncol(data_merged)])
colnames(data_merged) = c(names_files)
rownames(data_merged) = names_genes
write.table(data_merged, file="Counts_by_gene_total_for_DESeq2.txt", 
            row.names = T, sep="\t", quote=F)


setwd("../../../Scripts")
