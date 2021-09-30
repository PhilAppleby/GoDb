# Load the mongodb metadata markers collection
#
import time
import os, sys
import logging
from optparse import OptionParser
from dbsnpfile import Dbsnpfile # File access via tabix helper methods
from vcfrecord import VCFrecord

start_time = time.time()

def load_chrom_map(fh):
  chrom_map = {}
  for line in fh:
    data = line.strip().split()
    chrom_map[data[0]] = data[1]
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

  dbsnpfile = Dbsnpfile()
  dbsnpfile.set_tabix_file(options.dbsnpfile)
  for line in sys.stdin:
    count += 1
    line = line.strip()
    if (line.startswith('#')):
      print line
    else:
      vcfr = VCFrecord(line)
      posn = vcfr.get_posn_as_int()
      try:
        dbsnprecs = dbsnpfile.get_dbsnp_file_record(options.dbsnpfile, chrom_map[options.chrom], posn)
      except:
        print "Chromosome not found in map", options.chrom
        exit()
      if len(dbsnprecs) > 0:
        vcfr.set_varid(dbsnpfile.get_rsid(dbsnprecs[0]))
      else:
        logging.info("NOT FOUND for %s at %d" % (options.chrom, posn))
      print vcfr.get_record()

  return count
#
# execution flow starts here
#
parser = OptionParser()
parser.add_option("-d", "--dbsnpfile", dest="dbsnpfile",
   help="Path to file containing dbsnp rsids", metavar="STR")
parser.add_option("-c", "--chrom", dest="chrom",
   help="Numeric chromosome", metavar="STR")
parser.add_option("-m", "--chrommap", dest="chrommap",
   help="Map Numeric chromosome to chromosome string in NCBI format - eg: NC_000019.9", metavar="STR")
parser.add_option("-l", "--logfile", dest="logfile",
  help="logfile", metavar="STR")

(options, args) = parser.parse_args()
logging.basicConfig(filename=options.logfile,level=logging.DEBUG)

rec_count = main(options)
logging.info("END: %f seconds output=%d" % (time.time() - start_time, rec_count))
