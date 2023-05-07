#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(argparser))
suppressPackageStartupMessages(library(dplyr))

# this function is used to get information about biais in joining to coding/signal ends and resection around RSS sites
# You need to run the Tlxbed intersect script from your tlx files (with duplicates if needed) and the bed file of the regions extended on both sides of the RSS 
# All regions in the bed file should have the same length, use for ex mm10_Igk_40bpextend.bed, mm10_Igk_150bpextend.bed or mm10_IgkV_150bpextend.bed only for V fragments

# all junctions from within these bed files will be pooled together depending on their position compared to the RSS position (both the distance to RSS and the orientation towards coding or signal end of the RSS)
# output : you will get a txt file with the number of junctions in both coding orientation and signal orientation per position along the RSS 
# a simila file with normalized values is also generated 
# the stats file gives info about the total number of junctions, and about the % of junctions within some pre-determined intervals around the RSS 

# IMPORTANT NOTE ON RSS POSITION : the RSS position is located in between position -1 and 0 in the final files. 
# All positions <=-1 are junctions on the Coding side of th RSS
# All positions >= 0 are junctions on the Signal side of th RSS

# IMPORTANT NOTE ON QUANTIFICATION STATS1: quantification within 1bp of the RSS means the % of junctions located exactly on the nucleotide adjacent to RSS cutting position, on both directions
# Quantification within 5bp means the % of junctions located in a 10bp length intrval surrounding the RSS, 5bp on both sides

# QUANTIFICATION STATS2: to get the resection window within which at least x% of the junctions are

### Create a parser
p <- arg_parser("Resection around RSS and biais in joining")

### Add command line arguments
p <- add_argument(p, "meta", help = "metadata file for the library to analyse", type = "character")
p <- add_argument(p, "tlxbed_intersect_folder", help = "folder with results from tlsbed_intersect with the bedfile with extended regions on both parts of the RSS", type = "character")
p <- add_argument(p, "bedfile", help = "bedfile with extended regions on both sides of the RSS, used for the tlxbed intersect", type = "character")
p <- add_argument(p, "extension", help = "extension between your libraries name and the \"final _ext_intersect.tlx\" from tlxbed intersect, could be for exemple \"_redo\" or \"_JoinT\" ", type = "character")
p <- add_argument(p, "--which", help = "choose library to analyze", type = "character")
p <- add_argument(p, "--stats1", help = "steps for quantification: windows of resection", default = "2,10,40", type = "character")
p <- add_argument(p, "--stats2", help = "steps for quantification: % of junctions", default = "50, 90, 99", type = "character")
p <- add_argument(p, "--stats_only", flag = TRUE, help = "to be used to get new quantification files, if the script was run already and the output files lib_number.txt and lib_pct.txt are already there")

### Parse the command line arguments
argv <- parse_args(p)

### Prepare the work for each library
ResectionRSS <- function(meta, tlxbed_intersect_folder, bedfile, extension, ...){
  whichlib <- argv$which
  stats_only <- if_else(argv$stats_only, TRUE, FALSE)
  meta_all <- read.csv(meta, header=TRUE, sep="\t")
  lib_to_process <- seq_len(nrow(meta_all))
  bedfile_extendRSS <- read.csv(bedfile, sep="\t", header=FALSE)
  colnames(bedfile_extendRSS) <- c("chr","start","end","name", " ","strand")
  # calculate the extension on both sides of the RSS
  N <- (bedfile_extendRSS[1,]$end -  bedfile_extendRSS[1,]$start) /2
  
  # which libraries to process ?
  if (is.na(whichlib)){
    lib_to_process <- seq_len(nrow(meta_all))
  }
  if (!is.na(whichlib)){
    lib_to_process <- as.integer(unlist(strsplit(whichlib, ",")))
  }
  
  # write the first lines of the stats output files
  qtif1 <- argv$stats1
  qtif1 <- as.integer(unlist(strsplit(qtif1, ",")))
  qtif1_out <- paste("% of junctions located within ", qtif1, "bp of the RSS, i.e. in a ", 2*qtif1, "bp window", sep="")
  write(c("Library", "Total Number of junctions", "Number of junctions oriented towards Coding Sides of the RSS", "Number of junctions oriented towards Signal Sides of the RSS", "Ratio CE/SE", qtif1_out), file="Stats1.txt", ncolumns = (5 + length(qtif1)), sep = "\t")
  qtif2 <- argv$stats2
  qtif2 <- as.numeric(unlist(strsplit(qtif2, ",")))
  qtif2_out <- paste("Resection length for ", qtif2, "% of the junctions (bp)", sep="")
  write(c("Library", "Total Number of junctions", "Number of junctions oriented towards Coding Sides of the RSS", "Number of junctions oriented towards Signal Sides of the RSS", "Ratio CE/SE", "Resection Mean", "Resection Median", qtif2_out), file="Stats2.txt", ncolumns = (7 + length(qtif2)), sep = "\t")
  
  # for each library, collect the data
  for (i in lib_to_process) {
    lib <- paste(meta_all$Library[i], "_", meta_all$Sequencing[i], sep="")

    # collect the data from the tlx, unless --quantif_only is asked  
    if (stats_only){
      output_pct <- read.csv(paste(lib, "_pct.txt", sep =""), sep="\t", header =TRUE, row.names =1)
      output_nb <- read.csv(paste(lib, "_number.txt", sep =""), sep="\t", header =TRUE, row.names =1)
      nb_CE <- sum(output_nb[,1])
      nb_SE <- sum(output_nb[,2])
      tot_jcts <- sum(output_nb)
    } 
    
    else{
      print(paste("Collecting data for ", lib, sep=""))
      tlx <- paste(tlxbed_intersect_folder, "/", lib, extension, "_ext_intersect.tlx", sep="")
      data <- Position_collection(tlx, bedfile_extendRSS, N)
      tot_jcts <- sum(data$Number_of_junctions)
      output_nb <- matrix(data$Number_of_junctions, nrow = 2*N, ncol = 2, byrow=FALSE)
      rownames(output_nb) <- seq(-N, N-1, by=1)
      colnames(output_nb) <- c("Junctions oriented towards Coding sides of the RSS", "Junctions oriented towards Signal sides of the RSS")
      nb_CE <- sum(output_nb[,1])
      nb_SE <- sum(output_nb[,2])
      output_pct <- 100 * output_nb/tot_jcts
      write.table(output_nb, file=paste(lib, "_number.txt", sep =""), sep="\t", quote = FALSE, row.names = TRUE, col.names = NA)
      write.table(output_pct, file=paste(lib, "_pct.txt", sep =""), sep="\t", quote = FALSE, row.names = TRUE, col.names = NA)
    }
        
    # generate stats1 file : get % of junctions in defined intervals
    results1 <- c(lib, tot_jcts, nb_CE, nb_SE, round(nb_CE/nb_SE, digits = 4))
    for (x in qtif1){
      results1 <- c(results1, round(Quantif(output_pct, x, N), digits = 4))
    }
    write(results1, file="Stats1.txt", ncolumns = (5 + length(qtif1)), sep = "\t", append = TRUE)

    # generate stats2 file : get intervals for % of junctions
    output_df <- as.data.frame(output_pct)
    output_df$pos <- seq(-N, N-1, by=1)
    output_df <- output_df %>% mutate(res  = if_else( pos < 0, -pos - 1, pos ))

    # calculate the mean
    s <- 0
    for (k in (1:nrow(output_df))){
      s <- s + (output_df[k,1] + output_df[k,2]) * output_df$res[k] 
    }
    res_mean <- round(s / 100, digits = 4)
    # generate the output
    results2 <- c(lib, tot_jcts, nb_CE, nb_SE, round(nb_CE/nb_SE, digits = 4), res_mean, RSS_res(output_df, 50))
    for (x in qtif2){
      results2 <- c(results2, RSS_res(output_df, x))
    }
    write(results2, file="Stats2.txt", ncolumns = (7 + length(qtif2)), sep = "\t", append = TRUE)
  }
}


### Helpful functions

# function to get the beginning of the prey from one line of the tlx file, depending on the orientation
start_of_prey <- function(tlxfile, i){
  if (tlxfile$Strand[i] == -1) {
    return (tlxfile$Rend[i])
  }
  else {
    return (tlxfile$Rstart[i])
  }
}

# function to find the region from the bedfile where one junction from the tlxfile intersects
find_bed_line <- function(tlxfile, i, bedfile){
  for (l in 1:nrow(bedfile)){
    if (tlxfile$Rname[i] == bedfile$chr[l] && bedfile$start[l]-1 < start_of_prey(tlxfile,i) && bedfile$end[l]+1 > start_of_prey(tlxfile,i)){
      return(l)
    }
  }
}

# Function to quantification of the dispersion: 
# this gives the %of junction with n bp on both sides of the RSS, ie 2n bp around the RSS
# takes as argument the output_pct file with the % of junctions for each position
Quantif <- function(output_pct, n, N){
  if (n>N){
    print("Error: the range of quantification asked is larger than the extension used on both sides of the RSS")
    return()
  }
  sq <- as.character(seq(-n,n-1, by=1))
  return(sum(output_pct[sq,]))
}

# this gives the window of resection within which x% of the junctions are
RSS_res <- function(res_txt, pct){
  b <- 0
  j <- which(res_txt$pos >= -1 - b & res_txt$pos <= b)
  while ( (sum(res_txt[j,1] + res_txt[j,2])) < pct){
    b <- b + 1
    j <- which(res_txt$pos >= -1 - b & res_txt$pos <= b)
  }
  return(b)
}

# Function to collect positions of junctions for coding and signal ends in one single dataframe for each tlx
Position_collection <- function(tlx, bedfile_extendRSS, N){
  tlxfile <-read.csv(tlx, sep="\t")
  #create a data.frame to collect the positions of the junctions for signal ends and coding ends
  Ends_Collection = data.frame(Position = factor(1:(2*N)), Orientation = factor(c(rep("Coding",(2*N)),rep("Signal",(2*N)))), Number_of_junctions=rep(0,(2*N)))
  for (i in 1:nrow(tlxfile)){
    l <- find_bed_line(tlxfile, i, bedfile_extendRSS)
    BeginningOfPrey <- start_of_prey(tlxfile, i)
    #if the bedfile region for this junction is in positive orientation, it means that the position before RSS are the coding ends, the position after RSS are the signal ends
    #so if the junction is in positive orientation, it is a junction to the signal end
    #if the junction is in negative orientation, it is a junction to the coding end
    if (bedfile_extendRSS$strand[l] == "+") {
      position <- BeginningOfPrey-bedfile_extendRSS$start[l]
      if (tlxfile$Strand[i] == 1){
        Ends_Collection[which(Ends_Collection$Position == position & Ends_Collection$Orientation=="Signal"),]$Number_of_junctions <- 1 + Ends_Collection[which(Ends_Collection$Position == position & Ends_Collection$Orientation=="Signal"),]$Number_of_junctions
      }
      else {
        Ends_Collection[which(Ends_Collection$Position == position & Ends_Collection$Orientation=="Coding"),]$Number_of_junctions <- 1 + Ends_Collection[which(Ends_Collection$Position == position & Ends_Collection$Orientation=="Coding"),]$Number_of_junctions
      }
    }
    #if the bedfile region for this junction is in negative orientation, it means that the position before RSS are the signal ends, the position after RSS are the coding ends
    #so if the junction is in positive orientation, it is a junction to the coding end
    #if the junction is in negative orientation, it is a junction to the signal end
    if (bedfile_extendRSS$strand[l] == "-") {
      position <- bedfile_extendRSS$end[l] - BeginningOfPrey 
      if (tlxfile$Strand[i] == 1){
        Ends_Collection[which(Ends_Collection$Position==position & Ends_Collection$Orientation == "Coding"),]$Number_of_junctions <- 1 + Ends_Collection[which(Ends_Collection$Position == position & Ends_Collection$Orientation=="Coding"),]$Number_of_junctions
      }
      else {
        Ends_Collection[which(Ends_Collection$Position==position & Ends_Collection$Orientation == "Signal"),]$Number_of_junctions <- 1 + Ends_Collection[which(Ends_Collection$Position == position & Ends_Collection$Orientation=="Signal"),]$Number_of_junctions
      }
    }
  }
  return(Ends_Collection)
}

### Do work based on the passed arguments
ResectionRSS(argv$meta, argv$tlxbed_intersect_folder, argv$bedfile, argv$extension, argv$which, argv$stats1, argv$stats2, argv$stats_only)
