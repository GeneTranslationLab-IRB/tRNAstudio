
list.of.packages <- c("knitr","ggplot2","devtools","RColorBrewer", "ggrepel", "plyr","reshape2","FactoMineR",
                      "factoextra", "dplyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)){
  install.packages(new.packages, repos = "http://cran.us.r-project.org")
}

list.of.packages <- c("edgeR", "DESeq2","Rsamtools")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)>0){
  BiocManager::install(new.packages)
}

suppressPackageStartupMessages({
library("knitr")
library('Rsamtools')
library('ggplot2')
library('RColorBrewer')
library("DESeq2")
library("edgeR")
library('dplyr')
})

sample_data = read.table("..//Fastq_downloaded/sample_data.txt",
                         header=T, stringsAsFactors = T)
#sample_data = data.frame(ID=c("SRR7216347", "SRR7216348"), Condition=c("Control", "ADAT2"))
if (!file.exists("../Results/Plots")){
  dir.create("../Results/Plots", recursive=TRUE)
}

plot_data = c()



#1-Obtain count matrix.

tmp_total <- read.delim("../Results/R_files/Counts/Counts_by_gene_total_for_DESeq2.txt")
sum = apply(tmp_total,2,sum)
#2-Obtain sample info.

#3-Use the count matrix and the sample info to create the final matrix.
deseqdata_total <- DESeqDataSetFromMatrix(countData=tmp_total, 
                                          colData=sample_data, 
                                          design=~Condition)


## Sequencing depth

##Remember:The DESeq2 model internally corrects for library size, 
# so transformed or normalized values such as counts scaled by library size 
# should not be used as input.

#sum(colSums(assays(deseqdata_total)$counts)) #all the reads of the whole experiment 
#sum(colSums(assays(deseqdata_total)$counts))/ncol(deseqdata_total) #this is the ideal depth (is the same value as the mean of the summary)
dge <- DGEList(counts=assays(deseqdata_total)$counts, genes=mcols(deseqdata_total),
               group=deseqdata_total$Condition)

##Sequencing depth plot
#jpeg('../Results/Plots/Barplot_reads.jpeg', width=4100, height = 2900)
dge_plot = data.frame(Counts=dge$samples$lib.size/1e6, Group=dge$samples$group,
                      Sample =sample_data$ID )

p = ggplot(data=dge_plot, aes(y=Counts, x=Sample, fill=Group)) +
  geom_bar(stat="identity", width=0.75) + labs(y="Millions of counts", x="Samples") + 
  scale_fill_manual(values=c("darkred","navyblue"))+
  theme(axis.text=element_text(size=16, color="black"),panel.background = element_rect(fill = "grey92"),
      axis.title= element_text(size=20), legend.text = element_text(size = 16),
      legend.title = element_text(size=20), plot.title =element_text(size=20))

ggsave(p, file="Barplot_reads_total.jpeg", width = 35, 
       height = 25, units = "cm", path="../Results/Plots" )


for(condition in sample_data$Condition){
  dge_plot_group = dge_plot[dge_plot$Group == condition,]
  p = ggplot(data=dge_plot_group, aes(y=Counts, x=Sample, fill=Group)) +
    geom_bar(stat="identity", width=0.75) + labs(y="Millions of counts", x="Samples") +
    scale_fill_manual(values=c("darkred"))+
    theme(axis.text=element_text(size=18, color="black"),panel.background = element_rect(fill = "grey92"),
          axis.title= element_text(size=22), legend.text = element_text(size = 18),
          legend.title = element_text(size=22), plot.title =element_text(size=24))
  
  file_name = paste0("Barplot_reads_", condition, ".jpeg")
  ggsave(p, file=file_name, width = 35, 
         height = 25, units = "cm", path="../Results/Plots" )
}


## PCA plots (Samples clustering)
vsd <- varianceStabilizingTransformation(deseqdata_total, blind=FALSE)
rld <- rlog(deseqdata_total, blind=FALSE)
# this gives log2(n + 1)
ntd <- normTransform(deseqdata_total)
PCA_data = plotPCA(vsd, intgroup=c("Condition"), returnData= T)
plot = plotPCA(vsd, intgroup=c("Condition"))

label_x = plot$labels$x
label_y = plot$labels$y


p=ggplot(data=PCA_data, aes(x=PC1, y=PC2, fill=group)) +
  geom_point(size=10, shape=21) + 
  geom_text(label=PCA_data$name, vjust=3)+ 
  labs(x = label_x, y=label_y) + 
  guides(fill=guide_legend(title="Group")) + 
  scale_fill_manual(values=c("darkred","navyblue"))+
  theme(axis.text=element_text(size=15, color="black"),
        legend.text = element_text(size = 20), legend.title = element_text(size=20),
        axis.title= element_text(size=20))

file = "PCA_plot.jpeg"
ggsave(p, file=file, width = 35, 
       height = 25, units = "cm", path="../Results/Plots" )
#ggsave(PCA, file="PCA_plot2", width = 35,       height = 25, units = "cm", path="../Results/Plots" )

#Create an empty dataframe to marge all the sample resylts.
quality_data <- data.frame(Analysis= character(0), Counts= numeric(0), MAPQ= character(0))


for (sample in sample_data$ID) {
  #Precursor 
  precursor_mapped <- scanBam(paste0("../Results/",sample,"/Alignment_PG/",sample,"_PGloc_mapped.bam"))
  
  #Mature reads
  mature_mapped <- scanBam(paste0("../Results/",sample,"/Final_results/",sample,"all_mature.bam"))
  
  #Number of reads 
  precursor <- length(precursor_mapped[[1]]$qname)
  mature <- length(mature_mapped[[1]]$qname)
  total <- precursor+mature
  #MAPQ
  precursor_mapq <- as.data.frame(table(precursor_mapped[[1]]$mapq))
  mature_mapq <- as.data.frame(table(mature_mapped[[1]]$mapq))
  
  precursor_mapq$Var1 <-as.numeric(levels(precursor_mapq$Var1))[precursor_mapq$Var1]
  mature_mapq$Var1 <-as.numeric(levels(mature_mapq$Var1))[mature_mapq$Var1]
  
  #Select mapped reads according to the quality. 
  precursor_mapq_low <-subset(precursor_mapq,as.numeric(precursor_mapq$Var1) <= 2)
  precursor_mapq_good<-subset(precursor_mapq,as.numeric(precursor_mapq$Var1)  > 2)
  sum(precursor_mapq_low$Freq)+sum(precursor_mapq_good$Freq)
  
  mature_mapq_low <-subset(mature_mapq,as.numeric(mature_mapq$Var1) <= 2)
  mature_mapq_good <-subset(mature_mapq,as.numeric(mature_mapq$Var1)  > 2)
  sum(mature_mapq_low$Freq)+sum(mature_mapq_good$Freq)
  
  #Join precursor and mature results
  mapq_low <- sum(precursor_mapq_low$Freq)+sum(mature_mapq_low$Freq)
  mapq_good <- sum(precursor_mapq_good$Freq)+sum(mature_mapq_good$Freq)
  #Proportions
  quality_sample <- matrix(c(sample,mapq_low/total,'MAPQ 0-2',sample,mapq_good/total,'MAPQ > 2'),ncol=3,byrow=TRUE)
  #Add results to quality_data
  quality_data <- rbind(quality_data, quality_sample)
}

colnames(quality_data) <- c("Analysis","Counts","MAPQ")
quality_data$Counts<-as.numeric(quality_data$Counts)

######## Quality plot ########
quality_plot <- ggplot(quality_data, aes(fill=quality_data$MAPQ, y=quality_data$Counts*100, x=quality_data$Analysis)) +
  geom_bar(stat="identity",colour="black",width=0.65,lwd=0.3) + 
  theme_bw() +
  theme(legend.key.height = unit(0.5, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.text=element_text(size=12),
        text = element_text(size=12),
        axis.title.y = element_text(margin = margin(t = 0, r = 11, b = 0, l = 0)),
        axis.text.x = element_text(angle = 90,color="black", vjust = 0.5, hjust=1),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        plot.margin = unit(c(-0.75,0,0,0), "lines")) + 
  labs (x = "", y = "% of reads mapped to tRNAs") + 
  scale_fill_manual(values = c("#a9d491","grey90"), labels = c("MAPQ > 2","MAPQ <= 2")) +
  guides(fill=guide_legend(title="",size=5)) +  scale_y_reverse()

ggsave(quality_plot, file=paste0("../Results/Plots/Quality_barplot.jpeg"), 
       width = 20, height = 17, units = "cm")

for(group in sample_data$Condition){
  new_data = quality_data[quality_data$Analysis 
                          %in% sample_data$ID[sample_data$Condition ==group],]
  quality_plot <- ggplot(new_data, aes(fill=new_data$MAPQ, 
                                       y=new_data$Counts*100, x=new_data$Analysis)) +
    geom_bar(stat="identity",colour="black",width=0.65,lwd=0.3) + 
    theme_bw() +
    theme(legend.key.height = unit(0.5, "cm"),
          legend.key.width = unit(0.5, "cm"),
          legend.text=element_text(size=12),
          text = element_text(size=12),
          axis.title.y = element_text(margin = margin(t = 0, r = 11, b = 0, l = 0)),
          axis.text.x = element_text(angle = 90,color="black", vjust = 0.5, hjust=1),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          plot.margin = unit(c(-0.75,0,0,0), "lines")) + 
    labs (x = "", y = "% of reads mapped to tRNAs") + 
    scale_fill_manual(values = c("#a9d491","grey90"), labels = c("MAPQ > 2","MAPQ <= 2")) +
    guides(fill=guide_legend(title="",size=5)) +  scale_y_reverse()
  
  ggsave(quality_plot, file=paste0("../Results/Plots/Quality_barplot_",group, 
                            ".jpeg"),width = 20, height = 17, units = "cm")
}




