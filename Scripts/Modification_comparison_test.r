
list.of.packages <- c("GenomicAlignments","GenomicFeatures", "Rsamtools","edgeR",
                      "magick", "Glimma","reshape2", "RColorBrewer", "xlsx")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)){
  install.packages(new.packages, repos = "http://cran.us.r-project.org")
}


suppressPackageStartupMessages({
  library("Rsubread")
  library( "GenomicAlignments" )
  library( "GenomicFeatures" )
  library( "Rsamtools" )
  library("edgeR")
  library("vsn")
  library("magick")
  library("Glimma")
  library("pheatmap")
  library('tidyr')
  library('stringr')
  library(gtools)
  library(stringr)
  library(gdata)
  library(dplyr)
  library(plyr)
  library(data.table)
  library(heatmaply)
  library(gridExtra)
  library(ggplot2)
  library('RColorBrewer')
  library("xlsx")
})

join_log2fc_pvalue_data = function(data, log2fc, pvalue){
  log2fc = na.omit(log2fc)
  pvalue = na.omit(pvalue)
  if(length(log2fc) > 0){
    current_data = cbind(log2fc, pvalue)
    data = rbind(data, current_data)
  }
  data = data
  return(data)
}


sample_data = read.table("../Fastq_downloaded/sample_data.txt", header=T,
                         stringsAsFactors = T)

groups = levels(sample_data$Condition)
#groups = c("Circ", "DM")
dir = "../Results/Modification_ratio_plots"
dir_scripts = getwd()

aas = c("Ala", "Arg", "Asn", "Cys", "Gln", "Glu", "Gly", "His", "Ile", "iMet", 
        "Leu", "Lys", "Met", "Phe", "Pro", "SeC", "Ser", "Thr", "Trp", "Tyr", "Val")
files <- list.files(path = dir, pattern = ".txt")
if (!file.exists("../Results/Modification_ratio_plots/Comparison")){
  dir.create("../Results/Modification_ratio_plots/Comparison", recursive=TRUE)
}

aa_info = data.frame(pos="", gene="",pvalue="")
number_genes = 1
for(aa in aas){
  files_aa = grep(pattern=aa, files, value=T)
  genes_aa = levels(as.factor(grep(pattern=groups[1], files_aa, value=T)))
  genes_aa2 = levels(as.factor(grep(pattern=groups[2], files_aa, value=T)))
  
  gene_levels = substr(genes_aa, start=nchar(groups[1])+2, nchar(genes_aa))
  gene_levels2 = substr(genes_aa2, start=nchar(groups[2])+2, nchar(genes_aa2))
  # Only use the genes that have data in both groups.
  gene_levels = gene_levels[gene_levels %in% gene_levels2]
  
  if(aa %in% c("Ser", "Leu")){
    positions = c("9","26","32", "34", "37","45F", "58")
  }
  else{
    positions = c("9","26","32", "34", "37", "58")
  }
  position_info = matrix(ncol=length(positions), nrow=1)
  colnames(position_info) = positions
  rownames(position_info) = "pvalue"
  real_genes= substring(gene_levels, 1, nchar(gene_levels)-4)
  na =rep(NA,length(gene_levels))
  #
  nas = na
  for(i in 1:(length(positions)-1)){
    nas = cbind(nas,na)
  }
  gene_info = t(nas)
  colnames(gene_info) = real_genes
  rownames(gene_info) = positions
  cont = 0
  
  for(gene in real_genes){
    cont = cont +1
    data_group1 = read.table(file=paste0(dir,"/",groups[1],"_",gene, ".txt"), header=T)
    data_group2 = read.table(file=paste0(dir,"/",groups[2],"_",gene, ".txt"), header=T)
    colnames(data_group1)[2] = "Corrected_pos" 
    colnames(data_group2)[2] = "Corrected_pos" 
    
    data_group1 = data.frame(data_group1, Group = groups[1])
    data_group2 = data.frame(data_group2, Group =groups[2])
    data_tot = rbind(data_group1, data_group2)
    data_tot$Group = factor(data_tot$Group)
    p = ggplot(data = data_tot, aes(x=Corrected_pos, y=as.numeric(Modification_ratio), 
                                    group=Group, color = Group)) + geom_line() +
      labs(title='Modification ratio (Relative to gene coverage)', subtitle=gene,
           x="Position", y ="Modification ratio (Relative to gene coverage)")+
      scale_fill_manual(values= c("#0000FF", "#FF0000"),
                        breaks=c(groups[1], groups[2]))+
      scale_x_discrete(limits = as.character(data_group1$Corrected_pos))+
      theme(axis.text=element_text(size=16, color="black"),
            axis.text.x = element_text(size=12,angle=90),
            axis.title= element_text(size=20), legend.text = element_text(size = 16),
            legend.title = element_text(size=20), plot.subtitle = element_text(size=19),
            plot.title =element_text(size=20))
    
    file_name = paste0("Modification_ratio_by_pos_", gene, ".jpeg")
    ggsave(p, file=file_name, width = 35, 
           height = 25, units = "cm", path="../Results/Modification_ratio_plots/Comparison" )
    number_genes = number_genes+1
    ###FISHER TEST ####
    for(i in 1:length(positions)){
      if(positions[i] %in% data_group1$Corrected_pos){
        group1_pos<- as.data.frame(data_group1[data_group1$Corrected_pos == positions[i], ])
        group2_pos<- as.data.frame(data_group2[data_group2$Corrected_pos == positions[i], ])
        position_pos_test <- matrix(c(group1_pos$Base_coverage,group1_pos$Mod_bases,
                                      group2_pos$Base_coverage,group2_pos$Mod_bases),
                                    ncol=2,byrow=TRUE)
        
        rownames(position_pos_test)<-groups
        colnames(position_pos_test)<-c("Coverage","Modification")
        test_pos <- fisher.test(position_pos_test)
        p_pos <- test_pos$p.value
        position_info[i] = p_pos
      }
      else{
        position_info[i] = 1
      }
    }
    gene_info[,cont] = position_info
  }
  
  #without filtering the significative genes.
  #indx = gene_info<0.05
  #pos = rownames(gene_info)[row(gene_info)*indx]
  #gene = colnames(gene_info)[col(gene_info)*indx]
  #pvalue = gene_info[indx]
  pos = rownames(gene_info)
  gene = colnames(gene_info)
  pvalue = gene_info
  gene_info_final = data.frame(pos, gene, pvalue= as.numeric(pvalue))
  aa_info = rbind(aa_info, gene_info_final)
}


aa_info = aa_info[aa_info$pos!="",]
adjusted_pvalues = p.adjust(as.numeric(aa_info$pvalue), method="bonferroni")

aa_info_final = data.frame(aa_info, adjusted_pvalues)
  
#new_adjust = p.adjust(aa_info_final$pvalue[aa_info_final$pvalue<0.05], method = "bonferroni", 
                     # n = length(positions) * number_genes)
#aa_info_final$adjusted_pvalues[aa_info_final$pvalue<0.05] = new_adjust
#aa_info_final = aa_info_final[aa_info_final$adjusted_pvalues < 0.05,]


write.table(aa_info_final, file = "../Results/R_files/Modification_test.txt", 
            quote=FALSE, row.names = F)

#real_genes = substring(aa_info_final$gene, 1, nchar(aa_info_final$gene)-4)
unique_genes = unique(aa_info_final$gene)

pos_nine = c()
pos_twentysix = c()
pos_thirtytwo = c()
pos_thirtyfour = c()
pos_thirtyseven = c()
pos_fortyfiveF = c()
pos_fiftyeight = c()

for(aa in aas){
  id = grep(aa,unique_genes)
  if(aa %in% c("Ser", "Leu")){
    positions = c("9","26","32", "34", "37","45F", "58")
  }
  else{
    positions = c("9","26","32", "34", "37", "58")
  }
  genes_aa = unique_genes[id]
  final_data = matrix(ncol=length(positions), nrow=length(genes_aa))
  rownames(final_data) = genes_aa
  colnames(final_data) = positions
  
  pvalue_list = final_data
  mod1_list = final_data
  mod2_list = final_data
  cov1_list = final_data
  cov2_list = final_data
  custom_text = final_data
  relerror_list = final_data
  log2fc_list = final_data
  nt_mods_list1 = final_data
  nt_mods_list2 = final_data
  ref_list = final_data
  cont = 1 
  # Bucle for, in each isodecoder
  for(i in id){ 
    gene = unique_genes[i]
    data_group1 = read.table(file=paste0(dir,"/",groups[1],"_",gene, ".txt"),
                             header=T)
    data_group2 = read.table(file=paste0(dir,"/",groups[2],"_",gene, ".txt"), 
                             header=T)
    info_pvalue = c()
    info_mod1= c()
    info_mod2 = c()
    info_cov1 = c()
    info_cov2 = c()
    info_log2fc = c()
    info_nt_mods1 = c()
    info_nt_mods2 = c()
    info_ref = c()
    info_relerror = c()
    
    for(pos in positions){
      nas = TRUE
      if(pos %in% aa_info_final$pos[aa_info_final$gene == gene]){
        values1 = data_group1[data_group1$Consensus_tRNA_base_position ==pos,]
        values2 = data_group2[data_group2$Consensus_tRNA_base_position ==pos,]
        if(nrow(values1) >=1){
          cov1 = values1$Base_coverage
          cov2 = values2$Base_coverage
          if(cov1 > 50 & cov2 > 50){
            mod1 = values1$Modification_ratio * 100
            mod2 = values2$Modification_ratio * 100
            if(mod1 >10 | mod2>10){
              nas = FALSE
              
              rel_error = (mod2-mod1) / ((mod1+mod2)/2)
              log2fc = log2(mod2/mod1)
              
              ref = values1$Reference
              ref = substr(ref, 5,6)
              nt_mods1 = paste0("A: ", values1$A," G: ",  values1$G, " C: ", values1$C, 
                                " T: ", values1$T)
              nt_mods2 = paste0("A: ", values2$A," G: ",  values2$G, " C: ", values2$C, 
                                " T: ", values2$T)
              
              new_pvalue = (aa_info_final$adjusted_pvalues[aa_info_final$pos==pos & aa_info_final$gene==gene])
              new_pvalue = new_pvalue[1]  # Sometimes it returns the value repeated.
              new_pvalue = formatC(new_pvalue,format="e", digits = 2 )
              
              info_pvalue = c(info_pvalue, new_pvalue)
              info_mod1 = c(info_mod1, mod1)
              info_mod2 = c(info_mod2, mod2)
              info_cov1 = c(info_cov1, cov1)
              info_cov2 = c(info_cov2, cov2)
              info_log2fc = c(info_log2fc, log2fc)
              info_nt_mods1 = c(info_nt_mods1, nt_mods1)
              info_nt_mods2 = c(info_nt_mods2, nt_mods2)
              info_ref = c(info_ref, ref)
              info_relerror = c(info_relerror, rel_error)
            }
          }
        }
      }
      if(nas){
        info_pvalue = c(info_pvalue, NA)
        info_mod1 = c(info_mod1, NA)
        info_mod2 = c(info_mod2, NA)
        info_cov1 = c(info_cov1, NA)
        info_cov2 = c(info_cov2, NA)
        info_log2fc = c(info_log2fc, NA)
        info_relerror = c(info_relerror, NA)
        info_ref = c(info_ref, NA)
        info_nt_mods1 = c(info_nt_mods1, NA)
        info_nt_mods2 = c(info_nt_mods2, NA)
      }
    }
    
    final_data[cont,1:ncol(final_data)] = info_relerror
    mod1_list[cont,1:ncol(final_data)] = info_mod1
    mod2_list[cont,1:ncol(final_data)] = info_mod2
    pvalue_list[cont,1:ncol(final_data)] = info_pvalue
    cov1_list[cont,1:ncol(final_data)] = info_cov1
    cov2_list[cont,1:ncol(final_data)] = info_cov2 
    log2fc_list[cont,1:ncol(final_data)] = info_log2fc
    relerror_list[cont,1:ncol(final_data)] = info_relerror
    ref_list[cont,1:ncol(final_data)] = info_ref
    nt_mods_list1[cont,1:ncol(final_data)] = info_nt_mods1
    nt_mods_list2[cont,1:ncol(final_data)] = info_nt_mods2
    cont = cont + 1
  }
  
  if(nrow(final_data)>0){
    if(table(is.na(final_data)) != length(final_data)){
      setwd("../Results/Heatmaps")
      heatmap_file = paste0("Comparison_", aa, ".html")
      nas = apply(final_data, 1, function(x)  sum(is.na(x)))
      id = nas != ncol(final_data)
      names = colnames(final_data)
      final_data = round(data.frame(rbind(final_data[id,])), digits = 2)
      colnames(final_data) = names
      mod1_list = round(mod1_list[id,], digits = 2)
      mod2_list = round(mod2_list[id,], digits = 2)
      cov1_list = cov1_list[id,]
      cov2_list = cov2_list[id,]
      log2fc_list = round(log2fc_list[id,], digits = 2)
      pvalue_list = pvalue_list[id,]
      ref_list = ref_list[id,]
      nt_mods_list1 = nt_mods_list1[id,]
      nt_mods_list2 = nt_mods_list2[id,]
      relerror_list = round(relerror_list[id,], digits = 2)
      
      
      custom_text = data.frame(rbind(custom_text[id,]))
      custom_text[] = paste0("Gene: ", rownames(final_data), "\n",
                             "Modification ", groups[1], " (%): ", mod1_list, "\n",
                             "Modification ", groups[2], " (%): ", mod2_list, "\n",
                             "Coverage ", groups[1], ": ", cov1_list, "\n",
                             "Coverage ", groups[2], ": ", cov2_list, "\n",
                             "Reference ", ref_list, "\n",
                             groups[1], ": ", nt_mods_list1, "\n",
                             groups[2], ": ", nt_mods_list2, "\n",
                             "Relative difference: ", "(",groups[2], " - ", groups[1],")", 
                             " / mean = ", relerror_list, "\n", 
                             "Log2FC: ", log2fc_list, "\n",
                             "Adjusted pvalue: ", pvalue_list)
      
      #heatmap = heatmaply(final_data,  
      #                   colors= colorRampPalette(rev(brewer.pal(9, "RdBu"))),
      #                  plot_method = "plotly", 
      #                 limits=c(min(final_data, na.rm=TRUE),
      #                         max(final_data,na.rm=TRUE)),
      #               custom_hovertext=custom_text, Rowv = FALSE, Colv=FALSE, 
      #              xlab="Position", ylab="Gene", column_text_angle=0, 
      #             dendogram=FALSE, show_dendogram=c("FALSE", "FALSE"),
      #            file=heatmap_file)
      
      heatmap= heatmaply(final_data,
                         scale_fill_gradient_fun = 
                           ggplot2::scale_fill_gradient2(low = "blue",  high = "red", 
                                                         midpoint = 0,  limits= c(min(final_data, na.rm=TRUE),
                                                                                  max(final_data,na.rm=TRUE))),
                         custom_hovertext=custom_text, Rowv = FALSE, Colv=FALSE, 
                         xlab="Position", ylab="Gene", column_text_angle=0, 
                         dendogram=FALSE, show_dendogram=c("FALSE", "FALSE"),
                         file=heatmap_file)
      
      setwd(dir_scripts)
      

      if("matrix" %in% class(log2fc_list) ){
        pos_nine = 
          join_log2fc_pvalue_data(pos_nine, log2fc_list[,1], pvalue_list[,1])
        pos_twentysix = 
          join_log2fc_pvalue_data(pos_twentysix , log2fc_list[,2],pvalue_list[,2])
        pos_thirtytwo = 
          join_log2fc_pvalue_data(pos_thirtytwo, log2fc_list[,3],pvalue_list[,3])
        pos_thirtyfour = 
          join_log2fc_pvalue_data(pos_thirtyfour, log2fc_list[,4],pvalue_list[,4])
        pos_thirtyseven = 
          join_log2fc_pvalue_data(pos_thirtyseven, log2fc_list[,5],pvalue_list[,5])
        if("45F" %in% positions){
          pos_fortyfiveF = 
            join_log2fc_pvalue_data(pos_fortyfiveF, log2fc_list[,6],pvalue_list[,6])
          pos_fiftyeight = 
            join_log2fc_pvalue_data(pos_fiftyeight, log2fc_list[,7],pvalue_list[,7])
        }
        else{
          pos_fiftyeight = 
            join_log2fc_pvalue_data(pos_fiftyeight, log2fc_list[,6],pvalue_list[,6])
        }
      }  
      else{
        names_nine = rownames(pos_nine)
        pos_nine = 
          join_log2fc_pvalue_data(pos_nine, log2fc_list[1],pvalue_list[1])
        if(nrow(pos_nine)> length(names_nine)){
          rownames(pos_nine) = c(names_nine, gene)
        }
        names_twentysix = rownames(pos_twentysix)
        pos_twentysix = 
          join_log2fc_pvalue_data(pos_twentysix , log2fc_list[2],pvalue_list[2])
        if(nrow(pos_twentysix) > length(names_twentysix)){
          rownames(pos_twentysix) = c(names_twentysix, gene)
        }
        names_thirtytwo = rownames(pos_thirtytwo)
        pos_thirtytwo = 
          join_log2fc_pvalue_data(pos_thirtytwo, log2fc_list[3],pvalue_list[3])
        if(nrow(pos_thirtytwo) > length(names_thirtytwo)){
          rownames(pos_thirtytwo) = c(names_thirtytwo, gene)
        }
        names_thirtyfour = rownames(pos_thirtyfour)
        pos_thirtyfour = 
          join_log2fc_pvalue_data(pos_thirtyfour, log2fc_list[4],pvalue_list[4])
        if(nrow(pos_thirtyfour) > length(names_thirtyfour)){
          rownames(pos_thirtyfour) = c(names_thirtyfour, gene)
        }
        names_thirtyseven = rownames(pos_thirtyseven)
        pos_thirtyseven = 
          join_log2fc_pvalue_data(pos_thirtyseven, log2fc_list[5],pvalue_list[5])
        if(nrow(pos_thirtyseven) > length(names_thirtyseven)){
          rownames(pos_thirtyseven) = c(names_thirtyseven, gene)
        }
        if("45F" %in% positions){
          names_fortyfiveF = rownames(pos_fortyfiveF)
          pos_fortyfiveF = 
            join_log2fc_pvalue_data(pos_fortyfiveF, log2fc_list[6],pvalue_list[6])
          if(nrow(pos_fortyfiveF)>length(names_fortyfiveF)){
            rownames(pos_fortyfiveF) = c(names_fortyfiveF, gene)
          }
          names_fiftyeight = rownames(pos_fiftyeight)
          pos_fiftyeight = 
            join_log2fc_pvalue_data(pos_fiftyeight, log2fc_list[7],pvalue_list[7])
          if(nrow(pos_fiftyeight) > length(names_fiftyeight)){
            rownames(pos_fiftyeight) = c(names_fiftyeight, gene)
          }
        }
        else{
          names_fiftyeight = rownames(pos_fiftyeight)
          pos_fiftyeight = 
            join_log2fc_pvalue_data(pos_fiftyeight, log2fc_list[6],pvalue_list[6])
          if(nrow(pos_fiftyeight) > length(names_fiftyeight)){
            rownames(pos_fiftyeight) = c(names_fiftyeight, gene)
        }
      }  
      }
      
    }
  }
}


library("xlsx")
# Write the first data set in a new workbook

file_name = "../Results/R_files/Modification_comparison_results.xlsx"
write.xlsx(pos_nine, file = file_name,
           sheetName = "9", append = TRUE)
write.xlsx(pos_twentysix, file = file_name,
           sheetName = "26", append = TRUE)
write.xlsx(pos_thirtytwo, file = file_name,
           sheetName = "32", append = TRUE)
write.xlsx(pos_thirtyfour, file = file_name,
           sheetName = "34", append = TRUE)
write.xlsx(pos_thirtyseven, file = file_name,
           sheetName = "37", append = TRUE)
write.xlsx(pos_fortyfiveF, file = file_name,
           sheetName = "45F", append = TRUE)
write.xlsx(pos_fiftyeight, file = file_name,
           sheetName = "58", append = TRUE)
