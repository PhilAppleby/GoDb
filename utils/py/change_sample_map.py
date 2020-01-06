# Change sample conversions to a different code, omiyying to 
# change every nth one to force overlap checking
# 
import time
import os, sys
from optparse import OptionParser

start_time = time.time()

def main(options):
  intrval = int(options.number)

  count = 0
  for line in sys.stdin:
    line = line.strip()
    data = line.split(",")
    count += 1
    if count % intrval == 0:
      pass
    else:
      data[1] = options.prfx + data[1][1:]

    print ",".join(data)

  return count
#
# execution flow starts here
#
parser = OptionParser()
parser.add_option("-p", "--prfx", dest="prfx",
   help="Prefix for sample id's", metavar="STR")
parser.add_option("-n", "--number", dest="number",
   help="Interval for no-change to sample conversion", metavar="INT")

(options, args) = parser.parse_args()

if options.prfx == None:
  options.prfx = "2"
if options.number == None:
  options.number = "100"

rec_count = main(options)
#print "END:", time.time() - start_time, "seconds", rec_count
