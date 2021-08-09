
list.of.packages <- c("heatmaply","ggplot2", "RColorBrewer","dplyr","stringr" )
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)){
  install.packages(new.packages, repos = "http://cran.us.r-project.org")
}
suppressPackageStartupMessages({
library('heatmaply')
library('ggplot2')
library('RColorBrewer')
library('dplyr')
library('stringr')
})

order_positions = function(positions){
  loop = FALSE
  for(position in positions){
    if("V" %in% positions){
      loop =TRUE
    }
  }
  positions_ordered = positions[str_order(positions, numeric=TRUE)]
  if(loop){
    positions_loop = c("V11", "V12", "V13", "V14", "V15", "V16", "V17", "V1",
                       "V2", "V3","V4", "V5", "V27", "V26", "V25", "V24", "V23",
                       "V22", "V21")
    index_45 = match("45", positions_ordered)
    positions_ordered = append(positions_ordered, positions_loop,index_45)
  }
  positions_ordered
}


dir = "../Results/Modification_ratio_plots"

files <- list.files(path = dir, pattern = ".txt")

aas = c("Ala", "Arg", "Asn","Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile", "iMet", 
        "Leu", "Lys", "Met", "Phe", "Pro", "SeC", "Ser", "Thr", "Trp", "Tyr", "Val")

if (!file.exists("../Results/Heatmaps")){
  dir.create("../Results/Heatmaps", recursive=TRUE)
}

if (!file.exists("../Results/Reports")){
  dir.create("../Results/Reports", recursive=TRUE)
}

sample_data = read.table("..//Fastq_downloaded/sample_data.txt", header=T,
                         stringsAsFactors = T)

dir_scripts = getwd()
for(group in levels(sample_data$Condition)){
  files_group = grep(pattern=group, files, value=T)
  for(aa in aas){
    files_aa = grep(pattern=aa,files_group, value=T)
    
    data_modifications = data.frame()
    genes = c()
    for(l in 1:length(files_aa)){
      info = read.table(paste(dir,files_aa[l], sep="/"), header=TRUE)
      gene = substring(files_aa[l], (nchar(group)+2), nchar(files_aa[l])-4)
      
      genes = c(genes, gene)
      info = data.frame(gene, info)
      data_modifications = rbind(data_modifications, info)
    }
    
    ref <- data.frame(data_modifications,do.call(rbind,strsplit(as.character(data_modifications$gen),"-")))
    
    colnames(data_modifications)[3] = "Position_corrected"
    positions = unique(data_modifications$Position_corrected)
    
    positions_levels = order_positions(positions)
    data_modifications$Position_corrected = 
      factor(data_modifications$Position_corrected, levels=positions_levels)
    
    gene_levels = unique(data_modifications$gene)[
      str_order(unique(data_modifications$gene), numeric=TRUE)]
    data_modifications$gene = factor(data_modifications$gene, 
                                     levels=gene_levels)
    
    data_modifications_COVERAGE = data_modifications[data_modifications$Base_coverage > 50,]
    #levels(data_modifications$Position_corrected) = order(levels(data_modifications$Position_corrected)
    
    data_aa = matrix(ncol=length(levels(data_modifications$Position_corrected)),
                     nrow=length(levels(data_modifications$gene)))
    colnames(data_aa) = levels(data_modifications$Position_corrected)
    rownames(data_aa) = levels(data_modifications$gene)
    
    text_coverage = data_aa
    custom_text = data_aa
    reference = data_aa
    nt_mods = data_aa
    positions = data_aa
    
    data_aa_COVERAGE = data_aa
    custom_text_high_coverage = data_aa
    text_coverageCOVERAGE = data_aa
    reference_coverage = data_aa
    nt_mods_coverage = data_aa
    
    for(i in 1:ncol(data_aa)){
      for(j in 1:nrow(data_aa)){
        values = data_modifications %>% filter(gene==rownames(data_aa)[j] & Position_corrected== colnames(data_aa)[i]) %>% select(Modification_ratio, Base_coverage,Reference, Position_corrected)
        mods = data_modifications %>% filter(gene==rownames(data_aa)[j] & Position_corrected== colnames(data_aa)[i]) %>% select(A, G, C, T)
        
        values_coverage = data_modifications_COVERAGE %>% filter(gene==rownames(data_aa)[j] & Position_corrected== colnames(data_aa)[i]) %>% select(Modification_ratio, Base_coverage,Reference)
        mods_coverage = data_modifications_COVERAGE %>% filter(gene==rownames(data_aa)[j] & Position_corrected== colnames(data_aa)[i]) %>% select(A, G, C, T)
        positions[j,i] = colnames(data_aa)[i]
        
        #print(rownames(data_aa)[j])
        #print(colnames(data_aa)[i])
        
        if(nrow(values) !=0){
          
          data_aa[j,i] = values$Modification_ratio
          text_coverage[j,i] = values$Base_coverage
          reference[j,i] = substr(values$Reference, 5,6)
          
          names_mods = colnames(mods[which(mods!=0)])
          mods = data.frame(mods[,mods!=0])
          colnames(mods) = names_mods
          
          if(values$Base_coverage !=0){
            info_nt_mods = ""
            for(k in 1:ncol(mods)){
              info_nt_mods = paste(info_nt_mods, colnames(mods)[k], sep=" ")
              info_nt_mods = paste(info_nt_mods, mods[k], sep=": ")
            }
            info_nt_mods = substring(info_nt_mods, 2, nchar(info_nt_mods))
            nt_mods[j,i] = info_nt_mods
          }
        }
        
        if(nrow(values_coverage) !=0){
          data_aa_COVERAGE[j,i] = values_coverage$Modification_ratio
          text_coverageCOVERAGE[j,i] = values_coverage$Base_coverage
          reference_coverage[j,i] = values_coverage$Reference
          
          names_mods_coverage = colnames(mods_coverage[which(mods_coverage!=0)])
          mods_coverage = data.frame(mods_coverage[,mods_coverage!=0])
          colnames(mods_coverage) = names_mods_coverage
          
          if(values_coverage$Base_coverage !=0){
            info_nt_mods_coverage = ""
            for(k in 1:ncol(mods_coverage)){
              info_nt_mods_coverage = paste(info_nt_mods_coverage, colnames(mods_coverage)[k], sep=" ")
              info_nt_mods_coverage = paste(info_nt_mods_coverage, mods_coverage[k], sep=": ")
            }
            info_nt_mods_coverage = substring(info_nt_mods_coverage, 2, nchar(info_nt_mods_coverage))
            nt_mods_coverage[j,i] = info_nt_mods_coverage
          }
      }
      }
    }
    
    heatmap_file = paste0(group, "_heatmap_", aa,".html" )
    heatmap_file_coverage = paste0(group, "_heatmap_", aa,"_high_coverage.html" )

    custom_text[] = paste0("Family: ", rownames(data_aa), "\n",
      "Position: ", positions, "\n","Modification ratio: ",
                           round(data_aa, digits=2), "\n", "Coverage: ", text_coverage,
                  "\n", "Reference: ", reference, "\n", nt_mods)
    custom_text_high_coverage[] = paste0("Family: ", rownames(data_aa_COVERAGE), "\n",
      "Position: ", positions, "\n", "Modification ratio: ", 
      round(data_aa_COVERAGE, digits=2), "\n", "Coverage: ", text_coverageCOVERAGE,
                    "\n","Reference: ", reference_coverage, "\n", nt_mods_coverage)
    
    setwd("../Results/Heatmaps")
    
    heatmap = heatmaply(data_aa,  colors= colorRampPalette(brewer.pal(9, "OrRd")), plot_method = "plotly", limits=c(0,1),
                        Rowv = FALSE, Colv=FALSE, xlab="Position", ylab="Gene", custom_hovertext=custom_text, 
                        column_text_angle=0, dendogram=FALSE, show_dendogram=c("FALSE", "FALSE"), file=heatmap_file)
    
    
    if(table(is.na(data_aa_COVERAGE)) != dim(data_aa_COVERAGE)[1] * dim(data_aa_COVERAGE)[2]){
      heatmap = heatmaply(data_aa_COVERAGE,  colors= colorRampPalette(brewer.pal(9, "OrRd")), plot_method = "plotly", limits=c(0,1),
                        Rowv = FALSE, Colv=FALSE, xlab="Position", ylab="Gene", custom_hovertext=custom_text_high_coverage, 
                        column_text_angle=0, dendogram=FALSE, show_dendogram=c("FALSE", "FALSE"), file=heatmap_file_coverage)

      
      }
    
    #rm(heatmap)
    setwd(dir_scripts)
    #dev.off()
  }
  
  output_file = paste0("Report_", group)
  rmarkdown::render("Report.Rmd", params = list(group = group), 
                    output_file =output_file, output_dir = "../Results/Reports")
}