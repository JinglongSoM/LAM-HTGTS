# LAM-HTGTS
Containing Scripts for LAM-HTGTS data analysis
## Call TranslocPreprocess.pl in your terminal in the format as followings
TranslocPreprocess.pl path_to/metadata.txt path_to_preprocess/ --read1 path_to_sequence_file/Library_Run_Name.R1.fastq.gz --read2 path_to_sequence_file/Library_Run_Name.R2.fastq.gz 
### For advanced usage, TranslocPreprocess.pl --help

## Call TranslocPreprocess.pl in the following format
TranslocWrapper.pl path_to/metadata.txt preprocess/ results/ --threads N
### For advanced usage, TranslocWrapper.pl --help

## Call JoinT.R in the following format
path_to/JoinT.R path_to/meta.txt path_to/preprocess_folder/ path_to/result_folder/ path_to/output_directory/
### For advanced usage, JoinT.R --help
### the metadata should be constructed as for rejoin, i.e. with these columns:
Library	Sequencing	Researcher	Assembly	Chr	Start	End	Strand	Breakseq	Breaksite	MID	Primer	Adapter	Cutter	5nt_BaitEnd	5nt_PreyStart	Amplicon	Description

## Usage of tlx2bed.py as followings
tlx2bed.py [-h] -f TLXFILE [-o OUTPREFIX] [-t {bed,both,bedgraph}] [-g {mm9,hg19,hg38,mm10}] [--v3]

## The function of tlxbedintersect.py depends on two other scripts, tlx2BED.pl and pullTLXFromBED.pl.
Usage: python tlxbedintersect.py path_to/directory_with_tlx_files path_to/region_of_interest.bed

## Call tlx2BED-MACS.pl in the following format
tlx2BED-MACS.pl path_to/result_or_final_tlx Name_Of_New_File.bed

## Call JctStructure.R in the following format
JctStructure.R path_to/meta.txt path_to/folder_with_tlx_files extension_of_the_files nb_max_of_MH_insertion

## Call ResectionRSS.R in the following format
ResectionRSS.R path_to/meta.txt path_to/folder_with_results_from_tlxbed_intersect path_to/bedfile.bed extension_of_the_files
