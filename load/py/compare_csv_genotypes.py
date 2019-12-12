# Rewrite the a VCF headers line with subsitute sample_ids
# 
import time
import os, sys
from optparse import OptionParser

start_time = time.time()

def load_sample_map(fh):
  sample_map = {}
  for line in fh:
    data = line.strip().split(",")
    sample_map[data[1]] = data[0]

  return sample_map

def load_geno_map(fh):
  geno_map = {}
  hdr = fh.readline()
  print hdr.strip()
  for line in fh:
    data = line.strip().split(",")
    geno_map[data[0]] = data[1:]

  return geno_map

def main(options):
  try:
    fh = open(options.convfile, "r") 
    smap = load_sample_map(fh)
    fh2 = open(options.genofile, "r") 
    gmap = load_geno_map(fh2)
  except:
    print "Unexpected error:", sys.exc_info()[0]
    exit()

  hdr = []
  hdrlen = 0
  count = 0
  hdr = sys.stdin.readline()
  for line in sys.stdin:
    data = line.strip().split(",")
    if data[0] in smap:
      if smap[data[0]] in gmap:
        if gmap[smap[data[0]]] != data[1:]:
          print "mismatch %s,%s vs %s,%s " % (data[0], data[1:], smap[data[0]], gmap[smap[data[0]]])
        else:
          pass
          #print "match %s,%s vs %s,%s " % (data[0], data[1:], smap[data[0]], gmap[smap[data[0]]])
      else:
        print "ERROR 1"


  return count
#
# execution flow starts here
#
parser = OptionParser()
parser.add_option("-c", "--convfile", dest="convfile",
   help="sample_id conversion", metavar="STR")
parser.add_option("-g", "--genofile", dest="genofile",
   help="original genotype file", metavar="STR")

(options, args) = parser.parse_args()

rec_count = main(options)
#print "END:", time.time() - start_time, "seconds", rec_count
