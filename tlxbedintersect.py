#!/usr/bin/env python3
import os, sys

#inputs from teh command line
directory_name = sys.argv[1]
bed_file = sys.argv[2]

#goes through directory input and does full analysis (extract, extract complement, and count) for each tlx file in the directory
directory_list = os.listdir("%s" %(directory_name))
for item in directory_list:
	if ".tlx" in item and not(item.startswith(".")):
		tlx_file = directory_name+"/"+item
		tlx_bed = item.replace(".tlx",".bed")
		os.system("tlx2BED.pl %s %s" % (tlx_file, tlx_bed) )
		os.system("perl -pi -e 's/\r\n|\n|\r/\n/g' %s" % tlx_bed)

		#performs the extract option
		ext_output = item.replace(".tlx", "_ext_intersect.bed")
		ext_tlx_output = item.replace(".tlx", "_ext_intersect.tlx")
		os.system("bedtools intersect -u -a %s -b %s > %s" % (tlx_bed, bed_file, ext_output) )
		os.system("pullTLXFromBED.pl %s %s %s" % (tlx_file, ext_output, ext_tlx_output) )

		#performs the extract complement option
		comp_output = item.replace(".tlx", "_comp_intersect.bed")
		comp_tlx_output = item.replace(".tlx", "_comp_intersect.tlx")
		os.system("bedtools intersect -v -a %s -b %s > %s" % (tlx_bed, bed_file, comp_output) )
		os.system("pullTLXFromBED.pl %s %s %s" % (tlx_file, comp_output, comp_tlx_output) )

		#performs the count option
		count_output = item.replace(".tlx", "_count.bed")
		os.system("bedtools intersect -c -a %s -b %s > %s" % (bed_file, tlx_bed, count_output) )

		#removes intermediary files
		os.system("rm -f %s %s %s" % (tlx_bed, ext_output, comp_output) )

print("All done!")
