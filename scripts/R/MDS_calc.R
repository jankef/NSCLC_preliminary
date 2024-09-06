#
# Script name: MDS_calc.R
#
# Author: Florian Janke
# Last updated: 20240906
#
# Description:
#   (1) Calculates motif diversity scores from <extract_length-ends.py> output
#
# -------------------------------------------------------------------------------------------------------
library(stringr)
library(data.table)
library(optparse)
library(matrixStats)
library(data.table)
library(stringr)
library(sva)


#---------------
# S E T   P A R A M E T E R S
#---------------
option_list <- list(
  
  make_option(c("--input"), type = "character", default = NULL, help = "Path to directory containing the .tsv output files from <aberrant_positions.py> to be used."),
  make_option(c("--output"), type = "character", default = NULL, help = "Output file path (.tsv)."),
  make_option(c("--batch"), type = "character", default = NULL, help = ".tsv file with two columns: 01, Sample_IDs w/o extensions; 02, batch number")
  
)

opt <- parse_args(OptionParser(option_list=option_list))
options(scipen=0, stringsAsFactors=F)


#-----------
# F U N C T I O N S
#-----------
count_matrix <- function(files, gc_normalize, both_strands){
  
  #---
  # Function
  #---
  prop_table <- function(files, exclude){
    
    # Exclude sample from files according to <exclude>
    if(!is.null(exclude)) files <- files[!grepl(exclude, files)]
    
    # Exclude cases
    files <- files[grepl("CTRL", files)]
    
    list <- list()
    for(i in 1:length(files)){
      
      # Load sample
      table <- data.frame(fread(files[i]))
      table <- table[table$End == "left_seq", ]
      
      # Get probabilities for fragment length and GC content
      type <- c("Length", "GC")
      for(n in 1:length(type)){
        
        if(type[n] == "Length") tmp <- table %>% dplyr::group_by(Length) %>% dplyr::summarise(Count = sum(Count, na.rm = TRUE))
        if(type[n] == "GC") tmp <- table %>% dplyr::group_by(GC) %>% dplyr::summarise(Count = sum(Count, na.rm = TRUE))
        
        # Calculate probability
        tmp$Prop <- tmp$Count /sum(tmp$Count, na.rm = TRUE)
        tmp$Count <- NULL
        
        # Re-name column
        colnames(tmp)[2] <- paste0("Sample_", i)
        
        if(i == 1) list[[type[n]]] <- tmp
        if(i != 1) list[[type[n]]] <- merge(list[[type[n]]], tmp, by = type[n])
        
      }
    }
    
    # Calculate median across all controls
    list[["Length"]] <- data.frame(Length = list[["Length"]]$Length, Prop = rowMedians(as.matrix(list[["Length"]][, c(2:ncol(list[["Length"]]))]), na.rm = TRUE))
    list[["GC"]] <- data.frame(GC = list[["GC"]]$GC, Prop = rowMedians(as.matrix(list[["GC"]][, c(2:ncol(list[["GC"]]))]), na.rm = TRUE))
    
    # Return <list>
    return(list)
    
  }
  
  #---
  # Main script
  #---
  
  for(i in 1:length(files)){
    
    # Load file
    tmp <- data.frame(fread(files[i]))
    
    # Extract name of current sample
    if(i != 1) prev_type <- ifelse(str_sub(name, 1, 4) == "CTRL", "control", "case") else prev_type <- NULL
    name <- tools::file_path_sans_ext(basename(files[i]))
    
    # Analyze fragment-ends at both strands
    if(both_strands){
      
      tmp[tmp$End == "right_seq", ]$Base <- strReverse(chartr("ATGC","TACG", tmp[tmp$End == "right_seq", ]$Base))
      
    } else tmp <- tmp[tmp$End == "left_seq", ]
    
    # Sum right and left end
    if(gc_normalize){
      
      # Get probability table
      if(is.null(prev_type)) prop <- prop_table(files = files, exclude = name) else {
        
        if(prev_type == "control") prop <- prop_table(files = files, exclude = name)
        
      }
      
      # Remove low count rows
      tmp <- tmp[tmp$Count >=5, ]
      
      # Exclude rows not existing in <prop>
      tmp <- tmp[tmp$Length %in% prop$Length$Length, ]
      tmp <- tmp[tmp$GC %in% prop$GC$GC, ]
      
      # Merge <tmp> and <prop>
      tmp <- merge(tmp, prop$Length, by = "Length", )
      tmp <- merge(tmp, prop$GC, by = "GC")
      
      # Multiply probabilities
      tmp$Prop <- tmp[[6]] * tmp[[7]]
      tmp <- tmp[, -c(6:7)]
      
      # Weight counts
      tmp$Count_weighted <- tmp$Count /tmp$Prop
      
      # Summarize by <Motif>
      tmp <- tmp %>% dplyr::group_by(Motif) %>% dplyr::summarise(Count = sum(Count_weighted, na.rm = TRUE))
      
    } else tmp <- tmp %>% dplyr::group_by(Motif) %>% dplyr::summarise(Count = sum(Count, na.rm = TRUE))
    
    # Calculate proportion
    tmp$Count <- tmp$Count /sum(tmp$Count, na.rm = TRUE)
    
    # Add sample name
    colnames(tmp)[2] <- tools::file_path_sans_ext(basename(files[i]))
    
    # Append
    if(i == 1) matrix <- tmp
    if(i != 1) matrix <- merge(matrix, tmp, by = "Motif")
    if(i == length(files)){
      
      rownames(matrix) <- matrix$Motif
      matrix$Motif <- NULL
      
    }
    
    # Message
    message(i)
    
  }
  return(matrix)
}


#-----------
# M A I N   S C R I P T
#-----------
# 1) List files
files <- list.files(opt$input, full.names = TRUE, pattern = "\\.tsv$")

# 2) Create count matrix
table <- count_matrix(files = files, gc_normalize = FALSE, both_strands = FALSE)

# 3) Batch correct controls
batch <- data.frame(fread(opt$batch))
batch <- batch[order(colnames(table)), ]
table <- round(data.frame(prop.table(as.matrix(table), 2)) *30000000, 0)
table <- data.frame(ComBat_seq(as.matrix(table), batch = batch[[2]], group = NULL))
table <- data.frame(prop.table(as.matrix(table), 2))

# 4) Transform to long format
table$Base <- rownames(table)
table <- table[, c(ncol(table), 1:(ncol(table)-1))]
data <- reshape2::melt(table)
data$variable <- gsub("\\.", "-", data$variable)
colnames(data)[c(2:3)] <- c("Sample_ID", "Prop")

# 5) Calculate motif diversity score
sampID <- names(table(data$Sample_ID))
for(i in 1:length(sampID)){
  
  name <- sampID[i]
  
  patID <- substr(name, 1, nchar(name) - 6)
  
  group <- ifelse(str_sub(name, 1, 4) == "CTRL", "control", "case")
  
  tmp <- data[data$Sample_ID == sampID[i], ]
  mds.tmp <- -sum(tmp$Prop * log(tmp$Prop)) /log(nrow(tmp))
  tmp <- data.frame(Sample_ID = name, Patient_ID = patID, group = group, MDS = mds.tmp)
  
  if(i == 1) mds <- tmp
  if(i != 1) mds <- rbind(mds, tmp)
  
}

fwrite(mds, file = opt$output, sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)




