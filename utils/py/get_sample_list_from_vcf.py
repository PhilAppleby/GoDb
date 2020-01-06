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
  for line in sys.stdin:
    line = line.strip()
    if (line.startswith('##')):
      pass
    else:
      if (line.startswith('#')):
        vcfr = VCFrecord(line)
        prfx, sfx = vcfr.get_prfx_sfx()
        for samp in sfx:
          print samp
        break


  return count


# execution flow starts here
#
rec_count = main()
#print "END:", time.time() - start_time, "seconds", rec_count

