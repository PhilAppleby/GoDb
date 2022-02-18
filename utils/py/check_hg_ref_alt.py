# ----------------------------------------------------------------
# Command line args:
# ----------------------------------------------------------------
# 
import time
import logging
import pymongo
import re
import os, sys
from optparse import OptionParser

from hgdb import Hgdb # db helper methods

start_time = time.time()
output_data_size = 100
fa_rec_size = 50
genome_str_len = output_data_size * fa_rec_size

def main(options):
  print options.assaytype

  try:
    connection = pymongo.MongoClient("localhost")
  except IOError as e:
    print "I/O error({0}): {1}".format(e.errno, e.strerror)
    exit()
  except TypeError as e:
    print "Missing arguments ", e
    exit()
  except:
    print "Unexpected error:", sys.exc_info()[0]
    sys.exit()

  db = connection.hg
  hgdb = Hgdb(db, genome_str_len)
  rec_count=0

  for line in sys.stdin:
    data = line.strip().split(",")
    compare = True
    if options.assaytype != None:
      if data[1] != options.assaytype:
        compare = False
    if compare == True:
      rec_count += 1
      ltr = hgdb.get_single_position(int(data[2]), int(data[3])).upper()
      if ltr != data[4]:
        print "Err Varid = %s, plat=%s, chr=%s, posn = %s, Genome ltr = %s, ref = %s, alt = %s" % (data[0], data[1], data[2], data[3], ltr, data[4], data[5])


  return rec_count

# execution flow starts here
#
parser = OptionParser()

parser.add_option("-c", "--assaytype", dest="assaytype",
   help="assaytype", metavar="STR")
parser.add_option("-l", "--logfile", dest="logfile",
    help="logfile", metavar="STR")

(options, args) = parser.parse_args()
logging.basicConfig(filename=options.logfile,level=logging.DEBUG)
rec_count = main(options)
logging.info("END: %f seconds output=%d" % (time.time() - start_time, rec_count))

