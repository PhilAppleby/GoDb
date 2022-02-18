# ----------------------------------------------------------------
# Command line args:
# ----------------------------------------------------------------
# 
import time
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
  print(options.chromosome)
  print(options.posn)

  try:
    connection = pymongo.MongoClient("localhost")
  except IOError as e:
    #print "I/O error({0}): {1}".format(e.errno, e.strerror)
    exit()
  except TypeError as e:
    print("Missing arguments {0}", e)
    exit()
  except:
    print("Unexpected error: {0}", sys.exc_info()[0])
    sys.exit()

  db = connection.hg
  hgdb = Hgdb(db, genome_str_len)

  return hgdb.get_single_position(int(options.chromosome), int(options.posn))

# execution flow starts here
#
parser = OptionParser()

parser.add_option("-c", "--chromosome", dest="chromosome",
   help="chromosome", metavar="STR")
parser.add_option("-o", "--posn", dest="posn",
   help="BP posn in chromosome)", metavar="STR")

(options, args) = parser.parse_args()
ltr = main(options).upper()
print ltr  
#print "END:", time.time() - start_time, "seconds", ltr, genome_str_len
