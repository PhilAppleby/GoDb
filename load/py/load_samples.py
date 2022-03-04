# Load the mongodb metadata collection
#
import time
import os, sys
from optparse import OptionParser
from vcfrecord import VCFrecord
from godb import GoDb

start_time = time.time()
flush_at = 1000

def main(options):
  try:
    godb = GoDb()
  except:
    print "Unexpected error:", sys.exc_info()[0]
    exit()

  hdr = []
  count = 0
  for line in sys.stdin:
    line = line.strip()
    if (line.startswith('##')):
      pass
    else:
      if (line.startswith('#')):
        vcfr = VCFrecord(line)
        prf, sfx = vcfr.get_prfx_sfx()
        for idx, field in enumerate(sfx):
          count += 1
          godb.process_sample_detail(field, idx, options.assaytype)
          if (godb.get_samples_len() > flush_at):
            godb.flush_sample_buff()
        break

  godb.flush_sample_buff()
  print ""
  return count
#
# execution flow starts here
#
parser = OptionParser()
parser.add_option("-a", "--assaytype", dest="assaytype",
  help="Assay type (illumina, affy, etc)", metavar="STR")

(options, args) = parser.parse_args()

rec_count = main(options)
print("END: {0:.5f} seconds output={1}".format(time.time() - start_time, rec_count))
