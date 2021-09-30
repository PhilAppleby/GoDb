# Load the mongodb metadata markers collection
#
import time
import os, sys
import logging
from optparse import OptionParser
from vcfrecord import VCFrecord

start_time = time.time()

def load_chrom_map(fh):
  chrom_map = {}
  for line in fh:
    data = line.strip().split()
    chrom_map[data[1]] = data[0]
  return chrom_map

def main(options):
  hdr = []
  hdrlen = 0
  count = 0
  try:
    fh = open(options.chrommap)
    chrom_map = load_chrom_map(fh)
  except:
    print "Unable to open", options.chrommap
    exit()

  for line in sys.stdin:
    count += 1
    line = line.strip()
    if (line.startswith('#')):
      print line
    else:
      vcfr = VCFrecord(line)
      strchrom = vcfr.get_chr()
      try:
        vcfr.set_chr(chrom_map[strchrom])
      except:
        logging.info("Chromosome not found in map %s, %s" % (options.chrommap, strchrom))
        exit()
      print vcfr.get_record()

  return count
#
# execution flow starts here
#
parser = OptionParser()
parser.add_option("-m", "--chrommap", dest="chrommap",
   help="Map Numeric chromosome to chromosome string in NCBI format - eg: NC_000019.9", metavar="STR")
parser.add_option("-l", "--logfile", dest="logfile",
  help="logfile", metavar="STR")

(options, args) = parser.parse_args()
logging.basicConfig(filename=options.logfile,level=logging.DEBUG)

rec_count = main(options)
logging.info("END: %f seconds output=%d" % (time.time() - start_time, rec_count))
