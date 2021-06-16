#! /usr/bin/Rscript


list.of.packages <- c("devtools","ggrepel", "stringr","BiocManager","ggplot2","plyr","RColorBrewer","dplyr","knitr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)){
  install.packages(new.packages, repos = "http://cran.us.r-project.org")
}

list.of.packages <- c("GenomeInfoDbData","Rsamtools", "edgeR","vsn","Glimma","ReportingTools")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)>0){
  
  BiocManager::install(new.packages, update=TRUE, ask=FALSE)
}

suppressPackageStartupMessages({
library('knitr')
library('Rsamtools')
library('ggplot2')
library('devtools')
library('RColorBrewer')
library('ggrepel')
library('plyr')
library('devtools')
library('stringr')
library('dplyr')
library('stringr')
})

options(scipen=999)
options(warn=-1)

#Variable with the group of samples to be analyzed and the name of the group.
sample_data = read.table("../Fastq_downloaded/sample_data.txt", header=T,
                         stringsAsFactors = T)

groups = levels(sample_data$Condition)

if (!file.exists("../Results/Counts_plots/Total/By_tRNA_family")){
  dir.create("../Results/Counts_plots/Total/By_tRNA_family", recursive=TRUE)
}

#sample_data = data.frame(ID=c("SRR7216347", "SRR7216348"), Condition=c("Control", "ADAT2"))

isoaceptor_tot = data.frame(Aminoacid="", Counts="", Condition="",stringsAsFactors=FALSE)
isodecoder_tot = data.frame(Codon="",Counts="", Condition="",stringsAsFactors=FALSE)
for(group in groups){
  
  trna_total=c()
  for(sample in sample_data$ID){
    if(sample_data$Condition[sample_data$ID==sample] == group){
      trna_new = read.delim(paste0("../Results/R_files/Counts/",sample,"all_total.txt"),row.names = NULL,header=F)
      #trna_new = read.delim(paste0("../Scripts/Counts/",sample,"all_total.txt"),row.names = NULL,header=F)
      
      trna_total = bind_rows(trna_total, trna_new)
    }
  }
  
  #Marge results.
  trna_total_all <- data.frame(aggregate(cbind(trna_total$V2), by=list(trna_total$V1), FUN=sum))
  trna_total <- data.frame(trna_total_all,do.call(rbind,strsplit(as.character(trna_total_all$Group.1),"-")))
  trna_total$X2 = factor(trna_total$X2)
  trna_total = trna_total[substring(trna_total$Group.1 ,1,2) == "tR",]
  
  #Create plots for each family classified by amino acid 

  aa <- levels(trna_total$X2)
  aa <- aa[aa !="*"]
  
  for (e in aa) {
    trna_total_sub <- subset(trna_total,trna_total$X2 ==e)
    
    trna_total_sub = trna_total_sub[str_order(trna_total_sub$Group.1, numeric=TRUE),]
    p <- ggplot (trna_total_sub, aes (y=as.numeric(V1), x=Group.1)) + 
      geom_bar (stat="identity",width=0.5, fill ='palegreen4') +
      theme (text = element_text(size=10),
             axis.text.x = element_text(angle=90, color="black", vjust=0.5), 
             axis.text.y = element_text(color="black"), 
             panel.background =element_rect(fill = "snow2", colour = "snow2",
                                            size = 0.5, linetype = "solid")) +
      labs (x = "tRNA family", y = "Counts")  +
      labs (title = "Total sequencing read counts", subtitle = e)
    ggsave (plot=p, filename=paste0("../Results/Counts_plots/Total/By_tRNA_family/",
                                    group,"_",e,"_by_family.jpeg"), width = 20, height = 10, units = "cm")
  }
  
  
  #Plots by anticodon
  trna_total$trna <- paste(trna_total$X2,trna_total$X3)
  trna_total_by_anticodon<-data.frame(aggregate(cbind(trna_total$V1), by=list(trna_total$trna), FUN=sum))
  trna_total_by_anticodon <- subset(trna_total_by_anticodon,trna_total_by_anticodon$Group.1 !="* *")


  p <- ggplot(trna_total_by_anticodon, aes(y=V1, x=Group.1 )) + 
    geom_bar (stat="identity",width=0.5, fill ='palegreen4') + 
    theme (text = element_text(size=10),axis.text.x = element_text(angle=90,color="black",vjust=0.5),
           axis.text.y = element_text(color="black")) + 
    labs (x = "tRNA isodecoder", y = "Counts") + 
    theme (panel.background = element_rect(fill = "snow2",colour = "snow2",size = 0.5, linetype = "solid")) + 
    labs (title = "Total sequencing read counts", subtitle = "Isodecoder")
  ggsave(plot=p, filename=paste0("../Results/Counts_plots/Total/",group,"_by_Isodecoder.jpeg"), width = 20, height = 10, units = "cm")
  trna_total_by_anticodon = data.frame(trna_total_by_anticodon, Condition = group)
  colnames(trna_total_by_anticodon) = c("Codon", "Counts", "Condition")
  
  isodecoder_tot = rbind(isodecoder_tot, trna_total_by_anticodon)
  
  #Plots by aminoacid 
  trna_total_by_aa <- data.frame(aggregate(cbind(trna_total$V1), by=list(trna_total$X2), FUN=sum))
  trna_total_by_aa <- subset(trna_total_by_aa,trna_total_by_aa$Group.1 !="*")

  
  p <- ggplot(trna_total_by_aa, aes(y=V1, x=Group.1 )) + 
    geom_bar(stat="identity",width=0.5, fill ='palegreen4') + 
    theme(text = element_text(size=10),axis.text.x = element_text(angle=90,color="black",vjust=0.5),
          axis.text.y = element_text(color="black")) + 
    labs(x = "tRNA isoacceptor", y = "Counts") +
    theme(panel.background = element_rect(fill = "snow2",colour = "snow2",size = 0.5, linetype = "solid")) +
    labs(title = "Total sequencing read counts", subtitle = "Isoacceptor")
  ggsave(plot=p, filename=paste0("../Results/Counts_plots/Total/",group,"_by_Isoacceptor.jpeg"), width = 20, height = 10, units = "cm")
  
  trna_total_by_aa = data.frame(trna_total_by_aa, Condition=group)
  colnames(trna_total_by_aa) = c("Aminoacid", "Counts", "Condition")
  isoaceptor_tot = rbind(isoaceptor_tot, trna_total_by_aa)
}

isodecoder_tot = isodecoder_tot[isodecoder_tot$Codon!="",]
isodecoder_tot$Counts = as.numeric(isodecoder_tot$Counts)
isoaceptor_tot = isoaceptor_tot[isoaceptor_tot$Aminoacid!="",]
isoaceptor_tot$Counts = as.numeric(isoaceptor_tot$Counts)

p = ggplot(isodecoder_tot, aes(y=Counts, x=Codon, fill=Condition )) + 
  geom_bar (stat="identity",width=0.5, position=position_dodge2()) + 
  scale_fill_manual(values=c("sienna4","indianred2")) +
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle=90,size=6.5,vjust=0.5,color="black"),
        axis.text.y = element_text(color="black")) + 
  labs (x = "tRNA isodecoder", y = "Counts") + 
  theme (panel.background = element_rect(fill = "snow2",colour = "snow2",
                                         size = 0.5, linetype = "solid")) + 
  labs (title = "Total sequencing read counts", subtitle = "Isodecoder")

ggsave(plot=p, filename=paste0("../Results/Counts_plots/Total/Total_by_Isodecoder.jpeg"), width = 20, height = 10, units = "cm")

p = ggplot(isoaceptor_tot, aes(y=Counts, x=Aminoacid, fill=Condition )) + 
  geom_bar (stat="identity",width=0.5, position=position_dodge2()) + 
  scale_fill_manual(values=c("sienna4","indianred2")) +
  theme (text = element_text(size=10),axis.text.x = element_text(angle=90,size=6.5,vjust=0.5,color="black"),axis.text.y = element_text(color="black")) + 
  labs (x = "tRNA isodecoder", y = "Counts") + 
  theme (panel.background = element_rect(fill = "snow2",colour = "snow2",size = 0.5, linetype = "solid")) + 
  labs (title = "Total sequencing read counts", subtitle = "Isodecoder")

ggsave(plot=p, filename=paste0("../Results/Counts_plots/Total/Total_by_Isoacceptor.jpeg"), width = 20, height = 10, units = "cm")



