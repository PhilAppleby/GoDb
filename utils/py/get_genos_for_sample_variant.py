# ----------------------------------------------------------------
# Command line args:
# 1) rsid
# ----------------------------------------------------------------
#
import time
import logging
import pymongo
import re
import os, sys
from optparse import OptionParser

from godb import GoDb # db helper methods
from vcfrecord import VCFrecord

start_time = time.time()

def load_snpfile_data(fh):
  snps = []
  for line in fh:
    line=line.strip()
    snps.append(line)
  return snps

def main(options):
  included_assaytypes = {"affy":1, "illumina":1, "broad":1, "metabo":1, "exome":1}
  godb = GoDb()


# Data structures
  atype_list = []
  atype_posns = {}
  marker_list = []
  rsid_assaytypes = {}
  rsid_dict = {}
  rsid_prfx_dict = {}
  rsid_cr_dict = {}
  rsid_info_dict = {}
  count = 0

# Step 1 - get the list of entries for each rsid - one per assaytype

  vardocs = godb.get_multiple_variants(options.rsid)

  sampposns = godb.get_sample_posns(options.sampleid)

  for doc in vardocs:
    filepath = godb.get_filepath(doc["assaytype"], doc["chromosome"])
    rec = godb.get_variant_file_data(filepath, doc["chromosome"], doc["position"])
    vcfr = VCFrecord(rec)
    prfx, sfx = vcfr.get_prfx_sfx()
    if doc["assaytype"] in sampposns:
      count += 1
      print("{0},{1},{2},{3},{4},{5},{6}".format(options.rsid, doc["chromosome"], doc["position"], options.sampleid, doc["assaytype"], sampposns[doc["assaytype"]], sfx[sampposns[doc["assaytype"]]]))
      # print "%s,%s,%s,%s,%s,%d,%s" % (options.rsid, doc["chromosome"], doc["position"], options.sampleid, doc["assaytype"], sampposns[doc["assaytype"]], sfx[sampposns[doc["assaytype"]]])

  return count


# execution flow starts here
#
parser = OptionParser()
parser.add_option("-r", "--rsid", dest="rsid",
   help="dbSNP rsid", metavar="STR")
parser.add_option("-f", "--sampleid", dest="sampleid",
   help="sample id", metavar="STR")
parser.add_option("-p", "--prfx", dest="prfx",
   help="VCF file data path prefix", metavar="STR")

(options, args) = parser.parse_args()
rec_count = main(options)
print("END: {0:.5f} seconds output={1}".format(time.time() - start_time, rec_count))
