# LAM-HTGTS
Containing Scripts for LAM-HTGTS data analysis
# Call TranslocPreprocess.pl in your terminal in the format as followings
TranslocPreprocess.pl /path/to/your/meta/metadata.txt /path/to/preprocess/ --read1 /home/sequence_data/Library_Run_Name/Library_Run_Name.R1.fastq.gz --read2 /home/sequence_data/Library_Run_Name/Library_Run_Name.R2.fastq.gz 
## For advanced usage, TranslocPreprocess.pl --help

# Call TranslocPreprocess.pl in the following format.
TranslocWrapper.pl /path/to/your/meta/metadata.txt preprocess/ results/ --threads N
## For advanced usage, TranslocWrapper.pl --help
