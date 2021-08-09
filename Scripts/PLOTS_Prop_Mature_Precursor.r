#! /usr/bin/Rscript

list.of.packages <- c("devtools","ggrepel", "stringr","knitr","ggplot2","devtools","RColorBrewer","plyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)){
  install.packages(new.packages, repos = "http://cran.us.r-project.org")
}

list.of.packages <- c("Rsamtools")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)>0){
  BiocManager::install(new.packages)
}

suppressPackageStartupMessages({
library('knitr')
library('Rsamtools')
library('ggplot2')
library('RColorBrewer')
library('ggrepel')
library('plyr')
library('devtools')
library('stringr')
library('dplyr')
})

#Variable with the group of samples to be analysed and the name of the group.
sample_data = read.table("../Fastq_downloaded/sample_data.txt", header=T,
                         stringsAsFactors = T)

groups = levels(sample_data$Condition)

if (!file.exists("../Results/Counts_plots/Mature_vs_Precursor/By_tRNA_family")){
  dir.create("../Results/Counts_plots/Mature_vs_Precursor/By_tRNA_family", recursive=TRUE)
}

trna_tot_all_samples = data.frame()
trna_mature_all_samples = data.frame()

counts_isoaceptors = data.frame()
mature_isoaceptors = data.frame()
for(group in groups){
  trna_total=c()
  trna_mature=c()
  
  for(sample in sample_data$ID){
    if(sample_data$Condition[sample_data$ID==sample] == group){
      trna_new = read.delim(paste0("../Results/R_files/Counts/",sample,"all_total.txt"),
                            row.names = NULL,header=F)
      trna_total = bind_rows(trna_total, trna_new)
      
      trna_mature_new = read.delim(paste0("../Results/R_files/Counts/", sample,"all_mature_sort.txt"),row.names = NULL,header=F)
      trna_mature <- bind_rows(trna_mature,trna_mature_new) 
    }
  }
    
  trna_total_all <- data.frame(aggregate(cbind(trna_total$V2), by=list(trna_total$V1), FUN=sum))
  trna_mature_all <- data.frame(aggregate(cbind(trna_mature$V3), by=list(trna_mature$V1), FUN=sum))
  
  #Prec vs Mature ratio
  trna_total <- data.frame(trna_total_all,do.call(rbind,strsplit(as.character(trna_total_all$Group.1),"-")))
  trna_total <- subset(trna_total,trna_total$X2 !="*")
  trna_mature <- data.frame(trna_mature_all,do.call(rbind,strsplit(as.character(trna_mature_all$Group.1),"-")))
  trna_mature <- subset(trna_mature,trna_mature$X2 !="*")
  trna_total$X2 = factor(trna_total$X2)
  
  
  
  trna_total_sample = data.frame(trna_total, group = group)
  trna_tot_all_samples = rbind(trna_tot_all_samples, trna_total_sample)
  trna_total_isoaceptor <- data.frame(aggregate(cbind(trna_tot_all_samples$V1), 
                                                by=list(trna_tot_all_samples$X2), FUN=sum))
  trna_total_isoaceptor = data.frame(trna_total_isoaceptor, group = group)
  counts_isoaceptors = rbind(counts_isoaceptors, trna_total_isoaceptor )
  
  trna_mature_sample = data.frame(trna_mature, group = group)
  trna_mature_all_samples = rbind(trna_mature_all_samples, trna_mature_sample)
  trna_mature_isoaceptor <- data.frame(aggregate(cbind(trna_mature_sample$V1), 
                                                by=list(trna_mature_sample$X2), FUN=sum))
  trna_mature_isoaceptor = data.frame(trna_mature_isoaceptor, group = group)
  mature_isoaceptors = rbind(mature_isoaceptors, trna_mature_isoaceptor)
  
  aa <- levels(trna_total$X2)
  aa<- aa[aa !="*"]
  
  #For each family plot by aa.
  for (e in aa) {
    trna_total_sub <- subset(trna_total,trna_total$X2 ==e)
    trna_mature_sub <- subset(trna_mature,trna_mature$X2 ==e)
    
    trna_prec_mature_group <- transform(trna_total_sub, mature=trna_mature_sub$V1)
    trna_prec_mature_group <- trna_prec_mature_group[c("Group.1","V1","mature")]
    colnames(trna_prec_mature_group)<- c("TRNA.NAME","total", "mature")
    
    #Obtain proportions.
    trna_prec_mature_group <- transform(trna_prec_mature_group, precursor=trna_prec_mature_group$total-trna_prec_mature_group$mature)
    trna_prec_mature_group <- transform(trna_prec_mature_group, precursor.ratio=(trna_prec_mature_group$precursor/trna_prec_mature_group$total))
    trna_prec_mature_group <- transform(trna_prec_mature_group, mature.ratio=(trna_prec_mature_group$mature/trna_prec_mature_group$total))
    
    #Prepare the data to create the plots.
    trna_mature_by_group <- subset(trna_prec_mature_group, select=c(TRNA.NAME,mature.ratio))
    trna_mature_by_group <- transform(trna_mature_by_group, type="Processed")
    trna_prec_by_group <- subset(trna_prec_mature_group, select=c(TRNA.NAME,precursor.ratio))
    trna_prec_by_group  <- transform(trna_prec_by_group , type="Precursor")
    colnames(trna_prec_by_group) <- c("TRNA.NAME","counts", "Type")
    colnames(trna_mature_by_group) <- c("TRNA.NAME","counts", "Type")
    plot_data_by_group <- rbind(trna_prec_by_group, trna_mature_by_group)
    plot_data_by_group <- subset(plot_data_by_group, TRNA.NAME!="* *")
    plot_data_by_group$Type = factor(plot_data_by_group$Type, levels=c("Precursor", "Processed"))
    plot_data_by_group = plot_data_by_group[str_order(plot_data_by_group$TRNA.NAME, numeric = TRUE),]
    
    #Create and save the plots.
    p <- ggplot(plot_data_by_group, aes(fill=Type, y=counts, x=TRNA.NAME)) + 
      geom_bar (position="stack", stat="identity",width=0.5, colour="black") + 
      theme (text = element_text(size=10),
             axis.text.x = element_text(angle=90, color="black", vjust=0.5),
             axis.text.y = element_text(color="black"), 
             panel.background = element_rect(fill = "snow2",colour = "snow2",
                                             size = 0.5, linetype = "solid"),
             legend.position = "top")+ 
      scale_fill_manual("", values = c("Processed" = "#bddef0", "Precursor" = "#7180a7"))+
      ylim (0,1) + 
      labs (x = "tRNA family", y = "Proportion") +
      labs (title = "pre-tRNA vs processed tRNA", subtitle = e) #+ scale_fill_discrete(name = "", labels = c("pre-tRNA","processed tRNA"))
    ggsave (plot=p, filename=paste0("../Results/Counts_plots/Mature_vs_Precursor/By_tRNA_family/",group,"_",e,"_by_family.jpeg"), width = 20, height = 10, units = "cm")
  }
  

  
  #By anticodon
  trna_mature$trna <- paste(trna_mature$X2,trna_mature$X3)
  trna_mature<-data.frame(aggregate(cbind(trna_mature$V1), by=list(trna_mature$trna), FUN=sum))
  
  trna_total$trna <- paste(trna_total$X2,trna_total$X3)
  trna_total<-data.frame(aggregate(cbind(trna_total$V1), by=list(trna_total$trna), FUN=sum))
  trna_prec_mature <- transform(trna_mature, Total=trna_total$V1)
  
  #Obtain proportions.
  colnames(trna_prec_mature) <- c("TRNA.NAME","mature", "total")
  trna_prec_mature <- transform(trna_prec_mature, precursor=trna_prec_mature$total-trna_prec_mature$mature)
  trna_prec_mature <- transform(trna_prec_mature, precursor.ratio=(trna_prec_mature$precursor/trna_prec_mature$total))
  trna_prec_mature <- transform(trna_prec_mature, mature.ratio=(trna_prec_mature$mature/trna_prec_mature$total))
  
  trna_mature <- subset(trna_prec_mature, select=c(TRNA.NAME,mature.ratio))
  trna_mature <- transform(trna_mature, type="Processed")
  
  trna_prec <- subset(trna_prec_mature, select=c(TRNA.NAME,precursor.ratio))
  trna_prec <- transform(trna_prec, type="Precursor")
  colnames(trna_prec) <- c("TRNA.NAME","counts", "Type")
  colnames(trna_mature) <- c("TRNA.NAME","counts", "Type")
  plot_data <- rbind(trna_prec, trna_mature)
  plot_data <- subset(plot_data, TRNA.NAME!="Ile GAT")
  trna_prec_mature <- subset(trna_prec_mature, TRNA.NAME!="Ile GAT")
  plot_data$Type = factor(plot_data$Type, levels=c("Precursor", "Processed"))
  plot_data = plot_data[str_order(plot_data$TRNA.NAME, numeric = TRUE),]
  
  p <- ggplot(plot_data, aes(fill=Type, y=counts, x=TRNA.NAME)) + 
    geom_bar (position="stack", stat="identity",width=0.5,colour="black") + 
    #theme_bw() + 
    theme (axis.text.x = element_text(angle=90,color="black", size = 7, vjust=0.5),
           axis.text.y = element_text(color="black")) + 
    ylim (0,1) + 
    labs (x = "tRNA isodecoder", y = "Proportion") + 
    theme (panel.background = element_rect(fill = "snow2",colour = "snow2",
                size = 0.5, linetype = "solid"),
           legend.position = "top") +
    labs (title = "pre-tRNA vs processed tRNA", subtitle = "Isodecoder") +
    scale_fill_manual("", values = c("Processed" = "#bddef0", "Precursor" = "#7180a7"))
  
  ggsave(plot=p, filename=paste0("../Results/Counts_plots/Mature_vs_Precursor/",group,"_by_Isodecoder.jpeg"), width = 20, height = 10, units = "cm")
  
  
  #By amino acid.
  trna_prec_mature <- subset(trna_prec_mature, select=c(TRNA.NAME,mature,total,precursor))
  trna_prec_mature <- data.frame(trna_prec_mature,do.call(rbind,strsplit(as.character(trna_prec_mature$TRNA.NAME)," ")))
  
  trna_prec_mature<-data.frame(aggregate(cbind(trna_prec_mature$mature,trna_prec_mature$precursor,trna_prec_mature$total), by=list(trna_prec_mature$X1), FUN=sum))
  colnames(trna_prec_mature)<- c("TRNA.NAME","mature", "precursor","total")
  
  #Obtain proportions.
  trna_prec_mature <- transform(trna_prec_mature,precursor.ratio=(trna_prec_mature$precursor/trna_prec_mature$total))
  trna_prec_mature <- transform(trna_prec_mature, mature.ratio=(trna_prec_mature$mature/trna_prec_mature$total))
  
  trna_mature_by_aa <- subset(trna_prec_mature, select=c(TRNA.NAME,mature.ratio))
  trna_mature_by_aa <- transform(trna_mature_by_aa, type="Processed")
  colnames(trna_mature_by_aa)<- c("TRNA.NAME","counts", "Type")
  
  trna_precursor_by_aa <- subset(trna_prec_mature, select=c(TRNA.NAME,precursor.ratio))
  trna_precursor_by_aa <- transform(trna_precursor_by_aa, type="Precursor")
  colnames(trna_precursor_by_aa)<- c("TRNA.NAME","counts", "Type")
  
  plot_data_by_aa <- rbind(trna_precursor_by_aa,trna_mature_by_aa)
  plot_data_by_aa <- subset(plot_data_by_aa, TRNA.NAME!="* *")
  plot_data_by_aa$Type = factor(plot_data_by_aa$Type, levels=c("Precursor", "Processed"))
  
  p <- ggplot(plot_data_by_aa, aes(fill=Type, y=counts, x=TRNA.NAME)) + 
    geom_bar (position="stack", stat="identity",width=0.5, colour="black") + 
    theme(text = element_text(size=10),
          axis.text.x = element_text(angle=90,color="black", vjust=0.5),
          axis.text.y = element_text(color="black"),
          panel.background = element_rect(fill = "snow2",colour = "snow2",size = 0.5, linetype = "solid"),
          legend.position = "top") + 
    ylim (0,1) + 
    labs (x = "tRNA isoacceptor", y = "Proportion") + 
    labs (title = "pre-tRNA vs processed tRNA", subtitle = "Isoacceptor") + 
    scale_fill_manual("", values = c("Processed" = "#bddef0", "Precursor" = "#7180a7"))
  ggsave(plot=p, filename=paste0("../Results/Counts_plots/Mature_vs_Precursor/",group,"_by_Isoacceptor.jpeg"), width = 20, height = 10, units = "cm")
}

plot_data = data.frame(gene = trna_tot_all_samples$Group.1, ratio_mature = 
            trna_mature_all_samples$V1 / trna_tot_all_samples$V1,
            group = trna_tot_all_samples$group, aa = trna_tot_all_samples$X2 )

for(e in aa){
  plot_aa = plot_data[plot_data$aa ==e,]
  plot_aa$group = factor(plot_aa$group)
  
  p <- ggplot(plot_aa, aes(y=ratio_mature, x=gene, fill=group)) + 
    geom_bar(position = position_dodge(), 
             stat="identity",width=0.5,color="black") + 
    scale_fill_manual("",values = c("#bddef0", "#7180a7"))+
    theme (axis.text.x = element_text(angle=90,color="black", size = 7, vjust=0.5),
           axis.text.y = element_text(color="black")) + 
    ylim (0,1) + 

    labs (x = "tRNA gene", y = "Mature ratio") + 
    theme (panel.background = element_rect(fill = "snow2",colour = "snow2",
                                           size = 0.5, linetype = "solid"),
           legend.position = "top") +
    labs (title = "Mature ratio")
  ggsave (plot=p, filename=paste0(
    "../Results/Counts_plots/Mature_vs_Precursor/By_tRNA_family/Comparison_"
    ,e,"_by_family.jpeg"), width = 20, height = 10, units = "cm")
  
}

plot_data_isoaceptors = data.frame(aa= counts_isoaceptors$Group.1, 
                  mature_ratio = mature_isoaceptors$V1 / counts_isoaceptors$V1,
                  group = counts_isoaceptors$group)

p <- ggplot(plot_data_isoaceptors, aes(y=mature_ratio, x=aa, fill=group)) + 
  geom_bar(position = position_dodge(), 
           stat="identity",width=0.5,color="black") + 
  scale_fill_manual("",values = c("#bddef0", "#7180a7"))+
  theme (axis.text.x = element_text(angle=90,color="black", size = 10, vjust=0.5),
         axis.text.y = element_text(color="black")) + 
  ylim (0,1) + 
  
  labs (x = "tRNA Isoaceptors", y = "Mature ratio") + 
  theme (panel.background = element_rect(fill = "snow2",colour = "snow2",
                                         size = 0.5, linetype = "solid"),
         legend.position = "top") +
  labs (title = "Mature ratio")
file_name = "../Results/Counts_plots/Mature_vs_Precursor/Comparison_by_isoacceptor.jpeg"
ggsave (plot=p, filename=file_name, width = 20, height = 10, units = "cm")

plot_data_groups = data.frame(counts_isoaceptors$Group.1, counts_isoaceptors$V1, 
                              counts_isoaceptors$group, mature_isoaceptors$V1)
colnames(plot_data_groups) = c("Isoaceptor", "Counts", "group", "Mature")
plot_data_groups = data.frame(aggregate(cbind(plot_data_groups$Counts, 
              plot_data_groups$Mature), by= list(plot_data_groups$group), FUN=sum))
colnames(plot_data_groups) = c("Group", "Counts", "Mature")

plot_data_groups = data.frame(group=plot_data_groups$Group, mature_ratio=
                              plot_data_groups$Mature/plot_data_groups$Counts)
p <- ggplot(plot_data_groups, aes(y=mature_ratio, x=group, fill=group)) + 
  geom_bar(position = position_dodge(), 
           stat="identity",width=0.5,color="black") +
  scale_fill_manual("",values = c("#bddef0", "#7180a7"))+
  theme (axis.text.x = element_text(angle=90,color="black", size = 10, vjust=0.5),
         axis.text.y = element_text(color="black")) + 
  ylim (0,1) + 
  
  labs (x = "Condition", y = "Mature ratio") + 
  theme (panel.background = element_rect(fill = "snow2",colour = "snow2",
                                         size = 0.5, linetype = "solid"),
         legend.position="none") +
  labs (title = "Mature ratio")
p
file_name = "../Results/Counts_plots/Mature_vs_Precursor/Comparison_by_group.jpeg"
ggsave (plot=p, filename=file_name, width = 20, height = 10, units = "cm")
