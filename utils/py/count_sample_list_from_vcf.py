# Load the mongodb metadata samples collection
#
import time
import re
import os, sys
import random
import json
from optparse import OptionParser
from vcfrecord import VCFrecord

start_time = time.time()

def main():
  count = 0
  rcount = 0
  gcount = 0
  for line in sys.stdin:
    line = line.strip()
    if (line.startswith('##')):
      pass
    elif (line.startswith('#')):
      rcount += 1
      vcfr = VCFrecord(line)
      prfx, sfx = vcfr.get_prfx_sfx()
      for samp in sfx:
        #print samp
        count += 1
    elif (rcount == 1):
      vcfr = VCFrecord(line)
      prfx, sfx = vcfr.get_prfx_sfx()
      for geno in sfx:
        gcount += 1
      break

  return count, gcount


# execution flow starts here
#
scount, gcount = main()
print("END: {0:.5f} seconds sampcount={1}, genocount={2}".format(time.time() - start_time, scount, gcount))
