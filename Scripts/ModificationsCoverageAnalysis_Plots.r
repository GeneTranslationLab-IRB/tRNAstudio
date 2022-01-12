#! /usr/bin/Rscript

##Modification plots.

suppressPackageStartupMessages({
library('gtools')
library('gdata')
library('plyr')
library('data.table')
library('gridExtra')
library('Rsamtools')
library('ggplot2')
library('devtools')
library('stringr')
library('dplyr')
})

options(warn=-1)
### TOTAL COUNTS MATRIX


#Variable with the group of samples to be analysed and the name of the group.
sample_data = read.delim("../Fastq_downloaded/sample_data.txt", header=T,
                         stringsAsFactors = T)

groups = levels(factor(sample_data$Condition))


if (!file.exists("../Results/Modification_Coverage")){
  dir.create("../Results/Modification_Coverage", recursive=TRUE)
}



for(group in groups){
  base_call=c()
  
  for(sample in sample_data$ID){
    if(sample_data$Condition[sample_data$ID==sample] == group){
      base_call_new = read.delim(paste0("../Results/R_files/Base_calling/",
                                        sample,"_base_calls_total.txt"),
                        row.names = NULL,header=T)
      base_call <- bind_rows(base_call,base_call_new) 
    }
  }
    
  #Marge results.
  total_count<-data.frame(aggregate(cbind((base_call$A),base_call$C,base_call$G,
                          base_call$T,base_call$REF.COUNTS),
                          by=list(TRNA.POS=base_call$TRNA.POS), FUN=sum))
  
  #Split tRNA pos to obtain specific info and rename columns.
  result <- data.frame(total_count,do.call(rbind,strsplit(as.character(total_count$TRNA.POS),":")))
  colnames(result)<- c("TRNA.POS", "A","C","G","T","REF.COUNTS","tRNA_fam","POS","REF","POS.CON")
  result$tRNA_fam = factor(result$tRNA_fam)
  result<-result[!(result$POS.CON=="74" | result$POS.CON=="75" | result$POS.CON=="76"),]
  #Obtain all the families of trna that we are working with. 
  fam_id <- unique(result$tRNA_fam)
  fam_id <- as.character(levels(unique(result$tRNA_fam)))
  
  #Create a data frame to save al the results obtained in a_ref.
  all_data<-data.frame()
  
  
  for (trna in fam_id) {
    a <- result[result$tRNA_fam == trna, ]
    a<-a[order(as.numeric(as.character(a$POS))),]
    a$sum = rowSums(a[ , c(2:5)], na.rm = T)
    a_ref <- data.frame(a,do.call(rbind,strsplit(as.character(a$REF),"-")))
    a_ref$REF <- NULL
    a_ref$X1 <- NULL
    a_ref <- transform(a_ref, Mod.prop.base=1-(a_ref$REF.COUNTS/a_ref$sum))
    a_ref <- replace(a_ref, is.na(a_ref), 0)
    a_ref_info <- subset(a_ref,select=c(TRNA.POS,tRNA_fam,POS,POS.CON,Mod.prop.base,REF.COUNTS,sum,A,C,G,T))
    all_data<- rbind(all_data, a_ref_info)
  }
  
  #Now all_data contains the proportion of modificacion for all the tRNAs.
  all_data <- data.frame(all_data,do.call(rbind,strsplit(as.character(all_data$TRNA.POS),"-")))
  all_data<- subset(all_data, select=c(TRNA.POS,tRNA_fam,POS,POS.CON,sum,
                                       REF.COUNTS,Mod.prop.base,A,C,G,T, X2))
  
  all_data<- transform(all_data, Mod.bases=all_data$sum-all_data$REF.COUNTS)
  all_data<- subset(all_data, select=c(TRNA.POS,tRNA_fam,POS,POS.CON,sum,
                                       REF.COUNTS,Mod.prop.base,A,C,G,T,X2))
  
  #Save the name of all the tRNA fam.
  trna_fam <- str_sort(unique(all_data$tRNA_fam), numeric=TRUE)
  aas = str_sort(unique(all_data$X2), numeric=TRUE)
  
  #Plots for each family.
  for(aa in aas){
    by_aa = as.data.frame(all_data[all_data$X2 ==aa, ])
    max_coverage = max(by_aa$sum)
    plots_cov_aa = c()
    if (!file.exists("../Results/Modification_Coverage/Coverage_Plots")){
      dir.create("../Results/Modification_Coverage/Coverage_Plots", recursive=TRUE)
    }
    
    if (!file.exists("../Results/Modification_Coverage/Modification_Coverage_Plots")){
      dir.create("../Results/Modification_Coverage/Modification_Coverage_Plots", recursive=TRUE)
    }

    trna_fam <- str_sort(unique(by_aa$tRNA_fam), numeric=TRUE)
    i = 0
    for (fam in trna_fam) {
      i = i +1
      by_fam <- as.data.frame(by_aa[by_aa$tRNA_fam == fam, ])
      by_fam <- transform(by_fam, Mod.bases=by_fam$sum-by_fam$REF.COUNTS)
      length_trna <- length(by_fam$TRNA.POS)
      by_fam <- data.frame(by_fam,do.call(rbind,strsplit(as.character(by_fam$TRNA.POS),":")))
      by_fam$X3 = by_fam$X3[]
      by_fam_table <- subset(by_fam,select=c(POS,POS.CON,Mod.prop.base,
                                             sum,Mod.bases,X3,A,G,C,T))
      colnames(by_fam_table) <- c("Position","Consensus_tRNA_base_position","Modification_ratio","Base_coverage","Mod_bases","Reference","A","G","C","T")
      write.table(by_fam_table,file=paste0("../Results/Modification_Coverage/Modification_Coverage_Plots/",group,"_",fam,".txt"), sep = "\t", row.names = FALSE)
      df <- read.delim(paste0("../Results/Modification_Coverage/Modification_Coverage_Plots/",group,"_",fam,".txt"),header = TRUE)
      
      #Plot for the modification relative to base coverage.
      plot_mod_base <- ggplot(df,aes(x=as.numeric(Position),y=as.numeric(paste(Modification_ratio)))) +
        geom_line (size=0.5,color="#78c5dc") + 
        ylim (0,1) +
        scale_x_discrete (limits = as.factor(df$Consensus_tRNA_base_position)) +
        theme (axis.text.x = element_text(size=5,angle=90),axis.text.y = element_text(size=5)) +
        labs (x="Consensus tRNA base position", y ="Modification ratio") +
        theme (axis.title=element_text(size=5)) +
        labs (title = fam, subtitle = "Modification Ratio (Relative to base coverage)")
      #Plot for coverage.
      max_val <- max(df$Base_coverage, na.rm = TRUE)
      plot_cov <- ggplot(df,aes(x=as.numeric(Position),y=as.numeric(paste(Base_coverage)))) + 
        geom_line (size=0.5,color="#78c5dc") +
        ylim (0,max_val) +
        scale_x_discrete (limits = as.factor(df$Consensus_tRNA_base_position)) +
        theme (axis.text.x = element_text(size=5, angle=90),axis.text.y = element_text(size=5)) +
        labs (x="Consensus tRNA base position", y ="Base Counts") +
        theme (axis.title=element_text(size=5)) + 
        geom_area (alpha = 0.7, fill="#d4f0f0") +
        labs (subtitle = "Base Coverage")
      
      #Save plots
      a<- grid.arrange(plot_mod_base, plot_cov, nrow=2)
      ggsave(a, file=paste0("../Results/Modification_Coverage/Modification_Coverage_Plots/",group,"_",
                            fam,".jpeg"), width = 20, height = 10, units = "cm")
      rm(a)
      name_fam = substr(fam,6, nchar(fam))
      plot_cov2 <- ggplot(df,aes(x=as.numeric(Position),y=as.numeric(paste(Base_coverage)))) + 
        geom_line (size=0.5,color="darkslategray") +
        scale_x_discrete(limits=df$Consensus_tRNA_base_position, breaks=c("5","35","58", "73"))+
        ylim(0,max_val) + labs (y =name_fam, x="Position")+
        theme (axis.title.y = element_text(size=5),
               axis.title.x = element_text(size=4), 
               axis.text.y = element_text(hjust = 1,size=5), 
               axis.text.x = element_text(size=5)) + 
        geom_area (alpha = 1, fill="paleturquoise2")
      
      plots_cov_aa[[i]] = plot_cov2
      dev.off()
    }
    if(length(plots_cov_aa)<4){
      number_cols = length(plots_cov_aa)
    }
    if(length(plots_cov_aa)<20 & length(plots_cov_aa)>4){
      number_cols = 4
    }
    if(length(plots_cov_aa) %in% c(6, 9, 5)){
      number_cols = 3
    }
    if(length(plots_cov_aa)>20){
      number_cols = 5
    }
    if(length(plots_cov_aa)>35){
      number_cols = 6
    }

    
    coverage_aa = do.call(grid.arrange, args=c(plots_cov_aa, ncol=number_cols,top="Coverage"))
    ggsave(coverage_aa, file=paste0("../Results/Modification_Coverage/Coverage_Plots/",
                                    group,"_Coverage_",aa,".jpeg"), width = 20, 
           height = 10, units = "cm")
    
    dev.off()
  }
}

file.remove("Rplots.pdf")
