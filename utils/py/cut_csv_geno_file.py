#
# Cuts columns by name from a csv file
#
# Always outputs the id field as the first one, followed by the named columns
#
import time
import datetime
import os, sys
import logging
import csv
from optparse import OptionParser

start_time = time.time()
sys.stdout.flush()

def get_hdr_map(hdr):
  hdrmap = {}

  for i, colname in enumerate(hdr):
    hdrmap[colname] = i

  return hdrmap

def main(options):
  csvreader = None
  count=0

  try:
    csvfile =  open(options.csvfile, "r")
    csvreader = csv.reader(csvfile)
  except IOError as e:
    logging.fatal("I/O error({0}): {1}".format(e.errno, e.strerror))
    sys.exit()
  except TypeError as e:
    logging.fatal("Missing arguments ", e)
    sys.exit()
  except:
    logging.fatal("Unexpected error:", sys.exc_info())
    sys.exit()

  cols = []
  outhdr = []
  prfxlen=options.prfxlen

  try:
    hdr = csvreader.next()
    hdrmap = get_hdr_map(hdr)
    hdr = hdr[prfxlen:]
    # hdr has now been shortened by prfxlen, includes only the genotype columns
    # sorted to allow later file concatenation
    hdr.sort()
    outhdr.append("FID")
    outhdr.append("IID")
    colnames = options.colnames.split(',')
    # Process the header record to capture column indices and build the output
    # header record
    for i, col in enumerate(hdr):
      if col in colnames:
        outhdr.append(col)
        #print i, col
        cols.append(hdrmap[col])
    if options.writehdr == "Y":
      print ",".join(outhdr)
    #print cols

    for row in csvreader:
      outrec=[]
      outrec.append(row[0])
      outrec.append(row[0])
      count += 1
      # iterate over each row element
#      for i,elem in enumerate(row):
#        if i in cols:
#          outrec.append(elem)
      for idx in cols:
        outrec.append(row[idx])
      print ",".join(outrec)
  except:
    logging.fatal("Unexpected error (2):i %s" % (sys.exc_info()[0]))
    sys.exit()

  return count, len(outhdr)
#
# execution flow starts here
#
parser = OptionParser()
# csv file containing detail columns
parser.add_option("-c", "--csvfile", dest="csvfile",
  help="csv file containing genotype detail data", metavar="FILE")
# col names are comma separated - no complaint is made if a colname doesn't exist in the data
parser.add_option("-n", "--colnames", dest="colnames",
  help="column name", metavar="STR")
# prefix length, skipped and then filled with contents of col[0]
parser.add_option("-p", "--prfxlen", dest="prfxlen",
  help="prefix length", metavar="INT")
# write header - write header? Aids in concatenation
parser.add_option("-w", "--writehdr", dest="writehdr",
  help="write header", metavar="CHAR")
# logging
parser.add_option("-l", "--logfile", dest="logfile",
  help="logfile", metavar="STR")


(options, args) = parser.parse_args()
if options.prfxlen is None:
  options.prfxlen = 1
else:
  options.prfxlen = int(options.prfxlen)
if options.writehdr is None:
  options.writehdr = "Y"
logging.basicConfig(filename=options.logfile,level=logging.DEBUG)
#
rec_count, colcount = main(options)
logging.info("END: %f seconds, %d, %d" % (time.time() - start_time, rec_count, colcount))
