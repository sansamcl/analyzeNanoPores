#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(stringr)
library("Biostrings")

sample_name <- args[1]

blat_results <- read.table(paste0(args[1],"_blat_results.txt"))
fastas <- readDNAStringSet(paste0(args[1],".fa"))

crisprSite <- "GACGCTGCTAGACTACCAGT" %>% DNAStringSet %>% {names(.) <- "crispr";.}
pcrFull_2490 <- "CAGTGAGTGGTGCTGTTTCCctgaaggaagggactaagggacggtggcgcgggcccggaccggggccccggggcggcggcacggccgatatggcatgctgtcacaaagtaatgctgctgctggacaccgcgggcggcgccgcccgccacagccgggtccggcgggccgccctgcgcctcctcacctatctgagttgccgattcggcctggccagggtccactgggccttcaagttctttgactcgcagggggcgcggagccggccgtcccgcgtgtctgacttccgcgagctggggtcccgctcgtgggaggactttgaggaggagctggaggccaggctcgaggatcgcgcccacctgcccggcccggcgcccagggccacccacacgcacggcgccctgatggagacgctgctagactaccagtgggaccggcccgagatcacgtcgcccacgaagccgatcctgcggagcagcgggaggagactgctggacgtggagagcgaggccaaggaggccgaggccgcgctcgggggcttggtgaacgccgtcttcctcctggccccctgtccgcactcgcagagggagctgctgcagttcgtgtctgggtgcgaggcccaggcccagcgcctgccgcccacccctaagcaggtgatggagaagttgttgcccaagagagtccgggaagtcatggtcgcccgaaaaatcaccttctactgggtggataccaccgaatggtctaaggtaaggaaggttactgtcgtctcagatggcgtgcacggtgctttccttgctaaccaGAACTTGGCAGCGTCCTTAG" %>%
  DNAStringSet %>% {names(.) <- "2490";.}
pcrFull_2492 <- "TCTCTCGGACCACAGAGACCcgcccaccgtcactccctccaccaatcagactgaagttttctccactgcccctccaatcaccgccaaaaatggaagttgctctcgtaacttacgtaaggcaagcagccaataggagcaaaccttggaaactagctcccccaatagcagcagttcctggagagcgcgggagcctttcgaccagggtcccaaaggaaagcagtgagtggtgctgtttccctgaaggaagggactaagggacggtggcgcgggcccggaccggggccccggggcggcggcacggccgatatggcatgctgtcacaaagtaatgctgctgctggacaccgcgggcggcgccgcccgccacagccgggtccggcgggccgccctgcgcctcctcacctatctgagttgccgattcggcctggccagggtccactgggccttcaagttctttgactcgcagggggcgcggagccggccgtcccgcgtgtctgacttccgcgagctggggtcccgctcgtgggaggactttgaggaggagctggaggccaggctcgaggatcgcgcccacctgcccggcccggcgcccagggccacccacacgcacggcgccctgatggagacgctgctagactaccagtgggaccggcccgagatcacgtcgcccacgaagccgatcctgcggagcagcgggaggagactgctggacgtggagagcgaggccaaggaggccgaggccgcgctcgggggcttggtgaacgccgtcttcctcctggccccctgtccgcactcgcagagggagctgctgcagttcgtgtctgggtgcgaggcccaggcccagcgcctgccgcccacccctaagcaggtgatggagaagttgttgcccaagagagtccgggaagtcatggtcgcccgaaaaatcaccttctactgggtggataccaccgaatggtctaaggtaaggaaggtTACTGTCGTCTCAGATGGCG" %>%
  DNAStringSet %>% {names(.) <- "2492";.}
crisprSite <- c(crisprSite,pcrFull_2490,pcrFull_2492)


df <- lapply(1:nrow(blat_results),function(i){
  blat_results_1 <- blat_results[i,]
  block_sizes <- blat_results_1$V19 %>% str_split(.,",") %>% unlist %>% as.numeric %>% .[!is.na(.)]
  block_starts <- blat_results_1$V20 %>% str_split(.,",") %>% unlist %>% as.numeric %>% .[!is.na(.)]
  block_ends <- block_starts + block_sizes
  ranges=paste(block_starts,block_ends,sep="_")
  data.frame(blat_results_1$V14,block_starts,block_ends,ranges)
}) %>% do.call("rbind",.)

completeCRISPRSite <- df %>%
  .[which(.$block_starts < 43),] %>%
  .[which(.$block_ends > 66),] %>%
  .$blat_results_1.V14

df2 <- lapply(1:nrow(blat_results),function(i){
  blat_results_1 <- blat_results[i,]
  block_sizes <- blat_results_1$V19 %>% str_split(.,",") %>% unlist %>% as.numeric %>% .[!is.na(.)]
  block_starts <- blat_results_1$V20 %>% str_split(.,",") %>% unlist %>% as.numeric %>% .[!is.na(.)]
  block_ends <- block_starts + block_sizes
  gaps=paste(block_ends[-length(block_ends)],block_starts[-1],sep="_")
  if(length(gaps) > 0){
    data.frame(blat_results_1$V14,gaps)
  } else {
  }
}) %>% do.call("rbind",.)

df3 <- str_split(df2$gaps, "_") %>% unlist %>% as.numeric %>% rbind

gaps_after_crispr <- which(df3[,1] > 66)
gaps_before_crispr <- which(df3[,2] < 43)

df2 <- df2[-c(gaps_after_crispr,gaps_before_crispr),]

block_starts <- table(df$block_starts) %>% .[which(.>1)] %>% .[rev(order(.))]
block_ends <- table(df$block_ends) %>% .[which(.>1)] %>% .[rev(order(.))]
ranges <- table(df$ranges) %>% .[which(.>1)] %>% .[rev(order(.))]

gaps <- table(df2$gaps) %>% .[which(.>1)] %>% .[rev(order(.))]

gaps2 <- paste(gaps,names(gaps), sep="_")

directory_name <- paste(sample_name,"gap_fastas",sep="_")
dir.create(directory_name, showWarnings = FALSE)
setwd(directory_name)

for (i in 1:length(gaps)){
  gap <- gaps[i] %>% names %>% as.character
  gap2 <- gaps2[i]
  samples <- df2$blat_results_1.V14[which(df2$gaps == gap)]
  fastas[samples] %>%
    c(crisprSite,.) %>%
    writeXStringSet(., paste0(gap2,".fasta"))
}

fastas[completeCRISPRSite] %>%
  c(crisprSite,.) %>%
  writeXStringSet(., "wt.fasta")
