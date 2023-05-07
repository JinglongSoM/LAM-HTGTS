#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(argparser))
suppressPackageStartupMessages(library(dplyr))

### Create a parser
p <- arg_parser("Jct_structure")

### Add command line arguments
p <- add_argument(p, "meta", help = "metadata file for the library to analyse", type = "character")
p <- add_argument(p, "tlx_folder", help = "folder with results from tlsbed_intersect with the bedfile with extended regions on both parts of the RSS", type = "character")
p <- add_argument(p, "extension", help = "extension between your libraries name and the \"final _ext_intersect.tlx\" from tlxbed intersect, could be for exemple \"_redo\" or \"_JoinT\" ", type = "character")
p <- add_argument(p, "N", help="window to check for the max number of bases of MH and insertions", type="integer")
p <- add_argument(p, "--which", help = "choose library to analyze", type = "character")
p <- add_argument(p, "--region", help = "collect information about the junctions only in this region, in the form chrX:11111-22222, both extremities being included", type = "character")
p <- add_argument(p, "--except_region", help = "collect information about the junctions outside in this region, in the form chrX:11111-22222, both extremities being included in the region", type = "character")

### Parse the command line arguments
argv <- parse_args(p)

# Prepare the work for each library
Jct_structure <- function(meta, tlx_folder, extension, N, ...){
  whichlib <- argv$which
  region <- if_else(is.na(argv$region), FALSE, TRUE)
  except_region <- if_else(is.na(argv$except_region), FALSE, TRUE)
  if (region == TRUE & except_region == TRUE){
    print("you defined both --region and --except_region argument : the script does not deal with that. Please use a tlxbed intersect first to get your junctions of intersect")
    return()
  }
  
  # which region to use to filter the tlx files ?
  if (region == TRUE){
    region_window <- unlist(strsplit(argv$region, ":"))
    chr <- region_window[1]
    coordinates <- unlist(strsplit(region_window[2], "-"))
    st <- coordinates[1]
    en <- coordinates[2]
  }
  if (except_region == TRUE){
    out_window <- unlist(strsplit(argv$except_region, ":"))
    chr_out <- out_window[1]
    coordinates_out <- unlist(strsplit(out_window[2], "-"))
    st_out <- coordinates_out[1]
    en_out <- coordinates_out[2]
  }
  
  # which libraries to process ?
  meta_all <- read.csv(meta, header=TRUE, sep="\t")
  if (is.na(whichlib)){
    lib_to_process <- seq_len(nrow(meta_all))
  }
  if (!is.na(whichlib)){
    lib_to_process <- as.integer(unlist(strsplit(whichlib, ",")))
  }
  
  # write a stats file to get the number of junctions for each library
  col2 <- if_else(region == FALSE & except_region == FALSE, "in the whole tlx", if_else(region == TRUE, paste("in the region : ", argv$region, sep= ""), paste("outside of the region : ", argv$except_region, sep= "")))
  write(c("Library", "Total number of junctions in the tlx", paste("Total number of junctions ", col2, sep= ""), paste("Including junctions with MH and insertion numbers <= ", N, sep = "")), file="Stats.txt", ncolumns = 4, sep = "\t")
  
  # collect data in dataframes for both the number of junctions and the percentage of junctions
  df_nb <- data.frame(MH_ins = c(seq(-N,N, by=1)))
  df_pct <- data.frame(MH_ins = c(seq(-N,N, by=1)))
  names(df_nb)[1] <- "MH/insertion"
  names(df_pct)[1] <- "MH/insertion"
  
  # Do the job for each library to process:
  for (i in (1:length(lib_to_process))) {
    lib <- paste(meta_all$Library[lib_to_process[i]], "_", meta_all$Sequencing[lib_to_process[i]], sep="")
    tlx <- read.csv(paste(tlx_folder, "/", lib, extension, ".tlx", sep=""), header=TRUE, sep="\t")
    tlx_fil1 <- tlx
    if (region == TRUE){
      tlx_fil1 <- filter(tlx, Rname == chr & Junction >= st & Junction <= en)
    }
    if (except_region == TRUE){
      tlx_fil1 <- filter(tlx, !(Rname == chr_out & Junction >= st_out & Junction <= en_out))
    }
    tlx_fil2 <- filter(tlx_fil1, Qstart-B_Qend-1 >=-N & Qstart-B_Qend-1 <= N)
    
    # create the histogram to get the data of the distribution of jct structure
    histo <- hist(tlx_fil2$Qstart - tlx_fil2$B_Qend - 1, breaks=seq(-(N+1), N, by=1))
    df_nb <- cbind(df_nb, c(histo$counts))
    df_pct <- cbind(df_pct, round(c(histo$density * 100), digits=4))
    names(df_nb)[i+1] <- lib
    names(df_pct)[i+1] <- lib
    
    # update the stat file
    write(c(lib, nrow(tlx), nrow(tlx_fil1), nrow(tlx_fil2)), file="Stats.txt", ncolumns = 4, sep = "\t", append = TRUE)
  }
  # write the dataframes with the results
  write.table(df_nb, file="JctStructure_nb.txt", sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  write.table(df_pct, file="JctStructure_pct.txt", sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
}


### Do work based on the passed arguments
Jct_structure(argv$meta, argv$tlx_folder, argv$extension, argv$N, argv$which, argv$region)
