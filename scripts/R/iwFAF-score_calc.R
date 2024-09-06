#
# Script name: iwFAF-score_calc.R
#
# Author: Florian Janke
# Last updated: 20240906
#
# Description:
#   (1) Calculates iwFAF scores from <aberrant_positions.py> output
#   (2) Assums control samples are marked by file names starting with 'CTRL'
#
# -------------------------------------------------------------------------------------------------------
library(stringr)
library(data.table)
library(optparse)
library(matrixStats)


#---------------
# S E T   P A R A M E T E R S
#---------------
option_list <- list(
  
  make_option(c("--input"), type = "character", default = NULL, help = "Path to directory containing the .tsv output files from <aberrant_positions.py> to be used."),
  make_option(c("--output"), type = "character", default = NULL, help = "Output file path (.tsv).")
  
)

opt <- parse_args(OptionParser(option_list=option_list))
options(scipen=0, stringsAsFactors=F)


#-----------
# F U N C T I O N S
#-----------
iwFAF_pon <- function(input, pon){
  
  # Initiate <list>
  list <- list()
  
  # Get FAF per GC content and per length
  type <- c("GC", "Length")
  for(n in 1:length(type)){
    
    # Summarize counts
    for(i in 1:length(pon)){
      
      # Load data
      data <- data.frame(fread(file.path(input, paste0(pon[i], ".tsv"))))
      
      # Summarize by GC content or fragment-length
      if(type[n] == "GC"){
        
        tmp_ab <- data %>% dplyr::group_by(GC) %>% dplyr::summarise(Aberrant = sum(Aberrant, na.rm = TRUE))
        tmp_non <- data %>% dplyr::group_by(GC) %>% dplyr::summarise(Non_aberrant = sum(Non_aberrant, na.rm = TRUE))
        tmp <- merge(tmp_ab, tmp_non, by = "GC")
        
      }
      if(type[n] == "Length"){
        
        tmp_ab <- data %>% dplyr::group_by(Length) %>% dplyr::summarise(Aberrant = sum(Aberrant, na.rm = TRUE))
        tmp_non <- data %>% dplyr::group_by(Length) %>% dplyr::summarise(Non_aberrant = sum(Non_aberrant, na.rm = TRUE))
        tmp <- merge(tmp_ab, tmp_non, by = "Length")
        
      }
      
      tmp$total <- rowSums(tmp[, c(2:3)])
      tmp$prob <- tmp$Aberrant / tmp$total
      tmp$Aberrant <- tmp$Non_aberrant <- tmp$total <- NULL
      
      # Adjust column names
      colnames(tmp) <- c("type", pon[i])
      
      # Append
      if(i == 1) prob <- tmp
      if(i != 1) prob <- merge(prob, tmp, by = c("type"), all = TRUE)
      
    }
    
    # Get median probability
    prob <- data.frame(type = prob$type, Prob = rowMedians(as.matrix(prob[, c(2:ncol(prob))]), na.rm = TRUE))
    
    # Adjust column names
    colnames(prob)[1] <- type[n]
    
    # Add to <list>
    list[[type[n]]] <- prob
    
  }
  
  # Return <list>
  return(list)
  
}


#-----------
# M A I N   S C R I P T
#-----------
# Loop through samples in <input>
files <- list.files(opt$input, full.names = TRUE, pattern = "\\.tsv$")


for(i in 1:length(files)){
  
  # Define panel-of-normals
  pon <- files[str_sub(tools::file_path_sans_ext(basename(files)), 1, 4) == "CTRL"]
  pon <- pon[!grepl(files[i], pon)]
  pon <- tools::file_path_sans_ext(basename(pon))
  
  # Obtain the proportions of aberrant reads per GC content and per fragment length
  prop <- iwFAF_pon(input = opt$input, pon = pon)
  
  # Load data
  tmp <- data.frame(fread(files[i]))
  
  # Add GC- and size-based probability
  tmp <- merge(tmp, prop$GC, by = "GC", all.x = TRUE)
  tmp <- merge(tmp, prop$Length, by = "Length", all.x = TRUE)
  colnames(tmp)[5:6] <- c("FAF_GC", "FAF_size")
  
  # Exclude extremes
  tmp <- tmp[tmp$GC >=10 & tmp$GC <=90, ]
  tmp <- tmp[tmp$Length >=90, ]
  
  # Normalize aberrant fragments
  tmp$Aberrant_norm <- log2(1 /(tmp$FAF_GC *tmp$FAF_size)) *tmp$Aberrant
  
  # Normalize non-aberrant fragments
  tmp$Non_aberrant_norm <- log2(1 /((1 -tmp$FAF_GC) *(1 -tmp$FAF_size))) *tmp$Non_aberrant
  
  # Calculate iwFAF score
  score <- sum(tmp$Aberrant_norm) /(sum(tmp$Aberrant_norm) +sum(tmp$Non_aberrant_norm))
  
  # Summarize
  name <- tools::file_path_sans_ext(basename(files[i]))
  group <- ifelse(str_sub(name, 1, 4) == "CTRL", "control", "case")
  tmp <- data.frame(Sample_ID = name, Patient_ID = substr(name, 1, nchar(name) - 6), group = group, score = score)
  
  # Append
  if(i == 1) res <- tmp
  if(i != 1) res <- rbind(res, tmp)
  message(i)
  
}


fwrite(res, file = opt$output, sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
