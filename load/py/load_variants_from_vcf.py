# Load the mongodb metadata markers collection
#
import time
import os, sys
from optparse import OptionParser
from godb import GoDb

start_time = time.time()
flush_at = 10000

def main(options):
  try:
    godb = GoDb()
  except:
    print "Unexpected error:", sys.exc_info()[0]
    exit()

  hdr = []
  hdrlen = 0
  count = 0
  for line in sys.stdin:
    line = line.strip()
    if (line.startswith('#')):
      pass
    else:
      godb.process_variant_detail_vcf(line, options.assaytype)
      count += 1
      if (godb.get_variants_len() >= flush_at):
        godb.flush_variant_buff()
        print ".", time.time() - start_time, "seconds", count

  godb.flush_variant_buff()
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
