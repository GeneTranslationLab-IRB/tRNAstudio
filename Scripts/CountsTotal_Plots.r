#! /usr/bin/Rscript


suppressPackageStartupMessages({
  library('knitr')
  library('Rsamtools')
  library('ggplot2')
  library('RColorBrewer')
  library('stringr')
  library('dplyr')
  library('stringr')
  library('tidyverse')
})

options(scipen=999)
options(warn=-1)
#Variable with the group of samples to be analyzed and the name of the group.
sample_data = read.delim("../Fastq_downloaded/sample_data.txt", header=T)

groups = levels(factor(sample_data$Condition))

if (!file.exists("../Results/Counts_Plots/Total/By_tRNA_family")){
  dir.create("../Results/Counts_Plots/Total/By_tRNA_family", recursive=TRUE)
}


isoaceptor_tot = data.frame(Aminoacid="", MAPQ="", Counts="", Condition="",stringsAsFactors=FALSE)
isodecoder_tot = data.frame(Codon="",MAPQ="", Counts="", Condition="",stringsAsFactors=FALSE)
mitocondrial_tot = data.frame(tRNAmt="", MAPQ= "", Counts="", Condition="", 
                              stringsAsFactors = FALSE)
for(group in groups){
  trna_total=c()
  mitoc_tRNA_total = c()
  for(sample in sample_data$ID){
    if(sample_data$Condition[sample_data$ID==sample] == group){
      trna_new = read.delim(paste0("../Results/",sample, "/Counts/", sample,
        "_counts_total.txt"),row.names = NULL,header=TRUE)
      mitoc_tRNA_new = read.delim(paste0("../Results/", sample, "/Counts/",sample,
                                     "_counts_mitochondrial.txt"),row.names = NULL,header=TRUE)
      mitoc_tRNA_total = bind_rows(mitoc_tRNA_total, mitoc_tRNA_new)
      trna_total = bind_rows(trna_total, trna_new)
    }
  }
  # remove tRNAmt from totaltRNA
  trna_total = trna_total[!grepl("tRNAmt", trna_total$tRNA_ID),]

  #Marge results.
  trna_total_all <- data.frame(aggregate(cbind(trna_total$CountsMAPQ.2, trna_total$CountsMAPQ.2.1),
                                         by=list(trna_total$tRNA_ID), FUN=sum))
  trna_total <- data.frame(
    trna_total_all,do.call(rbind,strsplit(as.character(trna_total_all$Group.1),
                                          "-")))
  trna_total$X2 = factor(trna_total$X2)
  trna_total = trna_total[substring(trna_total$Group.1 ,1,2) == "tR",]
  trna_total = trna_total[,c(1:3, 5, 6)]
  colnames(trna_total) = c("Gene","Good_counts", "Bad_counts",
                           "Aminoacid", "Anticodon")
  good_counts = data.frame(trna_total[,c(1,2,4,5)], Type = "GoodMAPQ")
  bad_counts = data.frame(trna_total[,c(1,3,4,5)], Type = "BadMAPQ")
  colnames(good_counts)[2] = "Counts"
  colnames(bad_counts)[2] = "Counts"
  trna_counts = rbind(good_counts, bad_counts)
  
  tRNAmt <- data.frame(
    aggregate(cbind(mitoc_tRNA_total$CountsMAPQ.2, mitoc_tRNA_total$CountsMAPQ.2.1), 
              by=list(mitoc_tRNA_total$tRNA_ID), FUN=sum))
  tRNAmt <- data.frame(
    tRNAmt,do.call(rbind,strsplit(as.character(tRNAmt$Group.1),"-")))
  tRNAmt$X2 = factor(tRNAmt$X2)
  tRNAmt = tRNAmt[,c(1:3, 5, 6)]
  colnames(tRNAmt) = c("Gene", "Good_counts", "Bad_counts", 
                       "Aminoacid", "Anticodon")
  good_counts = data.frame(tRNAmt[,c(1,2,4,5)], Type = "GoodMAPQ")
  bad_counts = data.frame(tRNAmt[,c(1,3,4,5)], Type = "BadMAPQ")
  colnames(good_counts)[2] = "Counts"
  colnames(bad_counts)[2] = "Counts"
  tRNAmt_counts = rbind(good_counts, bad_counts)
  
  #Create plots for each family classified by amino acid 

  aa <- levels(trna_total$Aminoacid)
  aa <- aa[aa !="*"]
  
  for (e in aa) {
    counts_isoaceptor <- subset(trna_counts,trna_counts$Aminoacid ==e)
    
    counts_isoaceptor = counts_isoaceptor[str_order(counts_isoaceptor$Gene, numeric=TRUE),]
    p <- ggplot (counts_isoaceptor, aes(y=as.numeric(Counts), x=Gene, fill=Type)) + 
      geom_bar(stat="identity",width=0.5) +
      theme (text = element_text(size=10),
             axis.text.x = element_text(angle=90, color="black", vjust=0.5), 
             axis.text.y = element_text(color="black"), 
             panel.background =element_rect(fill = "snow2", colour = "snow2",
                                            size = 0.5, linetype = "solid")) +
      labs (x = "tRNA family", y = "Counts")  +
      scale_fill_manual(values = c("grey55", "palegreen4"), name = "",
                        labels = c("MAPQ <= 2", "MAPQ > 2")) + 
      labs (title = "Total sequencing read counts", subtitle = e)
    ggsave (plot=p, filename=paste0("../Results/Counts_Plots/Total/By_tRNA_family/",
                                    group,"_",e,"_by_family.jpeg"), width = 20, height = 10, units = "cm")
  }
  
  
  #Plots by anticodon
  trna_counts$trna = paste(trna_counts$Aminoacid, trna_counts$Anticodon)
  trna_total_by_anticodon = data.frame(
    aggregate(cbind(trna_counts$Counts), by=list(trna_counts$trna, 
                                                 trna_counts$Type), FUN=sum))
  colnames(trna_total_by_anticodon) = c("Anticodon", "MAPQ", "Counts")

  p <- ggplot(trna_total_by_anticodon, aes(y=Counts, x= Anticodon, fill=MAPQ )) + 
    geom_bar (stat="identity",width=0.5) + 
    theme (text = element_text(size=10),axis.text.x = element_text(angle=90,color="black",vjust=0.5),
           axis.text.y = element_text(color="black")) + 
    labs (x = "tRNA isodecoder", y = "Counts") + 
    scale_fill_manual(values = c("grey55", "palegreen4"), name = "",
                      labels = c("MAPQ <= 2", "MAPQ > 2")) + 
    theme (panel.background = element_rect(fill = "snow2",colour = "snow2",size = 0.5, linetype = "solid")) + 
    labs (title = "Total sequencing read counts", subtitle = "Isodecoder")
  ggsave(plot=p, filename=paste0("../Results/Counts_Plots/Total/",group,"_by_Isodecoder.jpeg"), width = 20, height = 10, units = "cm")
  trna_total_by_anticodon = data.frame(trna_total_by_anticodon, Condition = group)
  colnames(trna_total_by_anticodon) = c("Codon","MAPQ", "Counts", "Condition")

  isodecoder_tot = rbind(isodecoder_tot, trna_total_by_anticodon)
  
  #Plots by aminoacid 
  trna_total_by_aa <- data.frame(
    aggregate(cbind(trna_counts$Counts), by=list(trna_counts$Aminoacid, 
                                                 trna_counts$Type), FUN=sum))
  colnames(trna_total_by_aa) = c("Aminoacid", "MAPQ", "Counts")
  
  p <- ggplot(trna_total_by_aa, aes(y=Counts, x=Aminoacid, fill=MAPQ )) + 
    geom_bar(stat="identity",width=0.5) + 
    theme(text = element_text(size=10),axis.text.x = element_text(angle=90,color="black",vjust=0.5),
          axis.text.y = element_text(color="black")) + 
    scale_fill_manual(values = c("grey55", "palegreen4"), name = "",
                      labels = c("MAPQ <= 2", "MAPQ > 2")) +  
    labs(x = "tRNA isoacceptor", y = "Counts") +
    theme(panel.background = element_rect(fill = "snow2",colour = "snow2",size = 0.5, linetype = "solid")) +
    labs(title = "Total sequencing read counts", subtitle = "Isoacceptor")
  ggsave(plot=p, filename=paste0("../Results/Counts_Plots/Total/",group,"_by_Isoacceptor.jpeg"), width = 20, height = 10, units = "cm")
  
  trna_total_by_aa = data.frame(trna_total_by_aa, Condition=group)
  colnames(trna_total_by_aa) = c("Aminoacid", "MAPQ", "Counts", "Condition")
  isoaceptor_tot = rbind(isoaceptor_tot, trna_total_by_aa)
  
  
  ## mitoc tRNA.
  tRNAmt_counts$trna =  paste(tRNAmt_counts$Aminoacid, tRNAmt_counts$Anticodon)
  colnames(tRNAmt_counts) = c("Gene", "Counts", "Aminoacid", "Anticodon", 
                              "MAPQ", "trna")
  p <- ggplot(tRNAmt_counts, aes(y=Counts, x=Gene, fill=MAPQ )) + 
    geom_bar(stat="identity",width=0.5) + 
    theme(text = element_text(size=10),axis.text.x = element_text(angle=90,color="black",vjust=0.5),
          axis.text.y = element_text(color="black")) + 
    labs(x = "tRNA gene", y = "Counts") +
    scale_fill_manual(values = c("grey55", "palegreen4"), name = "",
                      labels = c("MAPQ <= 2", "MAPQ > 2")) +    
    theme(panel.background = element_rect(fill = "snow2",colour = "snow2",size = 0.5, linetype = "solid")) +
    labs(title = "mitochondrial tRNA genes")
  ggsave(plot=p, filename=paste0("../Results/Counts_Plots/Total/",group,
                                 "_mitochondrialtRNA.jpeg"), width = 20, height = 10, units = "cm")
  
  tRNAmt_by_group = data.frame(tRNAmt_counts, Condition=group)
  tRNAmt_condition =  data.frame(tRNAmt_counts[,c(1,5,2)], Condition = group)
  colnames(tRNAmt_condition) = c("tRNAmt", "MAPQ", "Counts", "Condition")
  mitocondrial_tot = rbind(mitocondrial_tot, tRNAmt_condition)
  
}

isodecoder_tot = isodecoder_tot[isodecoder_tot$Codon!="",]
isodecoder_tot$Counts = as.numeric(isodecoder_tot$Counts)
isoaceptor_tot = isoaceptor_tot[isoaceptor_tot$Aminoacid!="",]
isoaceptor_tot$Counts = as.numeric(isoaceptor_tot$Counts)

isodecoder_tot = data.frame(
  aggregate(isodecoder_tot$Counts, by=list(
    isodecoder_tot$Codon,isodecoder_tot$Condition), FUN=sum))
colnames(isodecoder_tot) = c("Codon", "Condition", "Counts")
p = ggplot(isodecoder_tot, aes(y=Counts, x=Codon, fill=Condition )) + 
  geom_bar(stat="identity", position=position_dodge2()) + 
  #scale_fill_manual(values=c("sienna4","indianred2")) +
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle=90,size=6.5,vjust=0.5,color="black"),
        axis.text.y = element_text(color="black")) + 
  labs (x = "tRNA isodecoder", y = "Counts") + 
  theme (panel.background = element_rect(fill = "snow2",colour = "snow2",
                                         size = 0.5, linetype = "solid")) + 
  labs (title = "Total sequencing read counts", subtitle = "Isodecoder")

ggsave(plot=p, filename=paste0("../Results/Counts_Plots/Total/Comparison_by_Isodecoder.jpeg"), width = 20, height = 10, units = "cm")


isoaceptor_tot = data.frame(
  aggregate(isoaceptor_tot$Counts, by=list(
    isoaceptor_tot$Aminoacid,isoaceptor_tot$Condition), FUN=sum))
colnames(isoaceptor_tot) = c("Aminoacid", "Condition", "Counts")
p = ggplot(isoaceptor_tot, aes(y=Counts, x=Aminoacid, fill=Condition )) + 
  geom_bar (stat="identity",width=0.5, position=position_dodge2()) + 
  #scale_fill_manual(values=c("sienna4","indianred2")) +
  theme (text = element_text(size=10),axis.text.x = element_text(angle=90,size=6.5,vjust=0.5,color="black"),axis.text.y = element_text(color="black")) + 
  labs (x = "tRNA isodecoder", y = "Counts") + 
  theme (panel.background = element_rect(fill = "snow2",colour = "snow2",size = 0.5, linetype = "solid")) + 
  labs (title = "Total sequencing read counts", subtitle = "Isodecoder")

ggsave(plot=p, filename=paste0("../Results/Counts_Plots/Total/Comparison_by_Isoacceptor.jpeg"), width = 20, height = 10, units = "cm")


### Plot tRNAmt genes
mitocondrial_tot = mitocondrial_tot[mitocondrial_tot$tRNAmt!="",]
mitocondrial_tot$Counts = as.numeric(mitocondrial_tot$Counts)
mitocondrial_tot = data.frame(
  aggregate(mitocondrial_tot$Counts, by=list(
    mitocondrial_tot$tRNAmt,mitocondrial_tot$Condition), FUN=sum))
colnames(mitocondrial_tot) = c("tRNAmt", "Condition", "Counts")
p = ggplot(mitocondrial_tot, aes(y=Counts, x=tRNAmt, fill=Condition )) + 
  geom_bar (stat="identity",width=0.5, position=position_dodge2()) + 
  #scale_fill_manual(values=c("sienna4","indianred2")) +
  theme (text = element_text(size=10),axis.text.x = element_text(angle=90,size=6.5,vjust=0.5,color="black"),axis.text.y = element_text(color="black")) + 
  labs (x = "tRNA isodecoder", y = "Counts") + 
  theme (panel.background = element_rect(fill = "snow2",colour = "snow2",size = 0.5, linetype = "solid")) + 
  labs (title = "Total sequencing read counts", subtitle = "mitocondrial tRNA genes")

ggsave(plot=p, filename=paste0("../Results/Counts_Plots/Total/Comparison_mitochondrialtRNAs.jpeg"), width = 20, height = 10, units = "cm")

