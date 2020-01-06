# Rewrite the a VCF headers line with subsitute sample_ids
# 
import time
import os, sys
from optparse import OptionParser
from vcfrecord import VCFrecord

start_time = time.time()

def load_sample_map(fh):
  sample_map = {}
  for line in fh:
    data = line.strip().split(",")
    sample_map[data[0]] = data[1]

  return sample_map

def main(options):
  try:
    fh = open(options.convfile, "r") 
    smap = load_sample_map(fh)
  except:
    print "Unexpected error:", sys.exc_info()[0]
    exit()

  hdr = []
  hdrlen = 0
  count = 0
  for line in sys.stdin:
    line = line.strip()
    if (line.startswith('##')):
      print line
    else:
      if (line.startswith('#')):
        vcfr = VCFrecord(line)
        prfx, sfx = vcfr.get_prfx_sfx()
        for idx, elem in enumerate(sfx):
          sfx[idx] = smap[elem]
        print "\t".join(prfx) + "\t" + "\t".join(sfx)
      else:
        print line

  return count
#
# execution flow starts here
#
parser = OptionParser()
parser.add_option("-s", "--convfile", dest="convfile",
   help="sample_id conversion", metavar="STR")

(options, args) = parser.parse_args()

rec_count = main(options)
#print "END:", time.time() - start_time, "seconds", rec_count
