# LAM-HTGTS
Containing Scripts for LAM-HTGTS data analysis
## Call TranslocPreprocess.pl in your terminal in the format as followings
TranslocPreprocess.pl /path/to/your/meta/metadata.txt /path/to/preprocess/ --read1 /home/sequence_data/Library_Run_Name/Library_Run_Name.R1.fastq.gz --read2 /home/sequence_data/Library_Run_Name/Library_Run_Name.R2.fastq.gz 
### For advanced usage, TranslocPreprocess.pl --help

## Call TranslocPreprocess.pl in the following format
TranslocWrapper.pl /path/to/your/meta/metadata.txt preprocess/ results/ --threads N
### For advanced usage, TranslocWrapper.pl --help

## Call JoinT.R in the following format
/path/to/JoinT.R path/to/meta.txt path/to/preprocess_folder/ path/to/result_folder/ path/to/output_directory/
### For advanced usage, JoinT.R --help
### the metadata should be constructed as for rejoin, i.e. with these columns:
Library	Sequencing	Researcher	Assembly	Chr	Start	End	Strand	Breakseq	Breaksite	MID	Primer	Adapter	Cutter	5nt_BaitEnd	5nt_PreyStart	Amplicon	Description

## Usage of tlx2bed.py as followings
tlx2bed.py [-h] -f TLXFILE [-o OUTPREFIX] [-t {bed,both,bedgraph}] [-g {mm9,hg19,hg38,mm10}] [--v3]

