# ----------------------------------------------------------------
# Command line args:
# ----------------------------------------------------------------
# 
import time
import pymongo
import re
import os, sys
from optparse import OptionParser

from genomefile import Genomefile # File access via tabix helper methods

start_time = time.time()
output_data_size = 100
fa_rec_size = 50
genome_str_len = output_data_size * fa_rec_size
vcfpattern = "chr%.1d.txt.gz"

def main(options):
  print options.dir
  print options.chromosome
  print options.offset

  filepath = options.dir + "/" + vcfpattern % (int(options.chromosome))

  print filepath

  tbxf = Genomefile(genome_str_len)
  tbxf.set_tabix_file(filepath)

  recs = tbxf.get_genome_file_records(int(options.chromosome), int(options.offset)) 

  print recs

  return 


# execution flow starts here
#
parser = OptionParser()

parser.add_option("-d", "--dir", dest="dir",
   help="directory containign vcf files", metavar="STR")
parser.add_option("-c", "--chromosome", dest="chromosome",
   help="chromosome", metavar="STR")
parser.add_option("-o", "--offset", dest="offset",
   help="BP posn in chromosome)", metavar="STR")

(options, args) = parser.parse_args()
main(options)
print "END:", time.time() - start_time, "seconds", genome_str_len

