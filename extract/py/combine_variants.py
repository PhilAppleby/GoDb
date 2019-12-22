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
from multibuffermerge import Multibuffermerge
from vcfrecord import VCFrecord
from mafhelper import Mafhelper
from hwehelper import Hwehelper

start_time = time.time()

def load_snpfile_data(fh):
  snps = []
  for line in fh:
    line=line.strip()
    snps.append(line)
  return snps

def main(options):
  #supported_assaytypes = {"bigtest":1, "affy":1, "illumina":1, "broad":1, "metabo":1, "exome":1}
  supported_assaytypes = {"affy":1, "illumina":1, "broad":1, "metabo":1, "exome":1}
  #supported_assaytypes = {"affy":1, "illumina":1, "broad":1, "exome":1}
  #supported_assaytypes = {"affy":1, "illumina":1}
  #supported_assaytypes = {"affy":1}
  #supported_assaytypes = {"bigtest":1}
  #supported_assaytypes = {"biggertest":1}
  rsids = []
  godb = GoDb()

  try:
    if options.snpfile != None:
      fh = open(options.snpfile, "r") 
      rsids = load_snpfile_data(fh)
    else:
      rsids = options.rsids.split(",")
  except IOError as e:
    print "I/O error({0}): {1}".format(e.errno, e.strerror)
    exit()
  except TypeError as e:
    print "Missing arguments ", e
    exit()
  except:
    logging.info("Unexpected error: %s", str(sys.exc_info()))
    sys.exit()

# Step 0 - initialise db connection and instanciate helper objects
  mafh = Mafhelper()
  hweh = Hwehelper()
# Data structures
  atype_list = []
  atype_posns = {}
  marker_list = []
  rsid_assaytypes = {}
  rsid_dict = {}
  rsid_prfx_dict = {}
  rsid_cr_dict = {}
  rsid_info_dict = {}
  hdr_pref = ["#CHROM",  "POS", "ID",  "REF", "ALT", "QUAL",  "FILTER",  "INFO",  "FORMAT"]

# Step 1 - get the list of entries for each rsid - one per assaytype

  for rsid in rsids:
    #logging.info("Processing rsid = %s", rsid)
    docs = godb.get_multiple_variants(rsid)
    if docs.count() > 0:
      rsid_assaytypes[rsid] = []
    else:
      logging.info("RSID %s NOTFOUND", rsid)
  #print docs

  # Step 1a - collect assaytypes and marker documents
  # At this point we're establishing a list order which must be observed throughout.
    for doc in docs:
      #logging.info("%s", str(doc))
      if doc["assaytype"] not in supported_assaytypes:
        continue
      if doc["assaytype"] not in atype_list:
        atype_list.append(doc["assaytype"])
      rsid_assaytypes[rsid].append(doc)
  logging.info(str(atype_list))
# Step 2 - collect lists of prochis (sample ids) by assaytype
  prochi_list = [[]] * len(atype_list)
  for i, atype in enumerate(atype_list):
    atype_posns[atype] = i
    prochi_list[i] = godb.get_samples(atype)
    #logging.info("SAMP %d, %s, %s", i, atype, str(prochi_list[i]))

  mm = Multibuffermerge(prochi_list)

# Step 3 - get combined col_header positions
# combo is a dict {posn:colname}
  combo = mm.get_combined_positions()
  #print len(combo)
# combocol is a list [colname1, colname2, ..., colname] again we keep the order of this intact
  combocol = mm.get_combined_columns()
  
# Step 4 - for each variant by rsid
  for rsid in rsid_assaytypes:
    if rsid not in rsid_dict:
      rsid_prfx_dict[rsid] = [[]] * len(atype_list)
      rsid_dict[rsid] = [[]] * len(atype_list)
      rsid_cr_dict[rsid] = [[]] * len(atype_list)
      rsid_info_dict[rsid] = [[]] * len(atype_list)
    #print len(rsid_assaytypes[rsid])
    for doc in rsid_assaytypes[rsid]:
      if options.prfx != None:
        fpath = godb.get_full_filepath(doc["assaytype"], doc["chromosome"], options.prfx)
      else:
        fpath = godb.get_filepath(doc["assaytype"], doc["chromosome"])
      logging.info("Assaytype=%s fpath=%s", doc["assaytype"], fpath)

      result = godb.get_variant_file_data(fpath, doc["chromosome"], doc["position"])
      if result != None:
        vcfr = VCFrecord(result)
        varid = vcfr.get_varid()
        if varid == rsid:
          rec = result
          maf, ma, cr = mafh.get_maf_and_cr(vcfr)
          # TODO - ALSO check maf, also apply QC filter at individual record level
          rsid_cr_dict[doc["rsid"]][atype_posns[doc["assaytype"]]] = cr
          rsid_dict[doc["rsid"]][atype_posns[doc["assaytype"]]] = rec
          logging.info("%s (%s), maf=%s, ma=%s, cr=%s" % (doc["rsid"], doc["assaytype"], maf, ma, cr))
  
  #print combocol
# Step 5 - execute the merge process
  print "\t".join(hdr_pref + combocol)
  count = 0
  concordant = True
  for rsid in rsid_dict:
    if len(rsid_dict[rsid][0]) > 0:
      if options.check == 'Y':
        concordant = mm.check_concordancies(rsid_dict[rsid], atype_list, options.chipval)

      if concordant == True:
        comborec = mm.get_combined_array(rsid_dict[rsid], rsid_cr_dict[rsid], atype_list)
        vcfr = VCFrecord(rsid_dict[rsid][0])
        prfx,sfx = vcfr.get_prfx_sfx()
        if len(prfx) > 0:
          logging.info("PRFX = %s, for %s", str(prfx), rsid)
          prfx[8] += ":AT"
          outrec = prfx + comborec
          print "\t".join(outrec)
          count += 1
        else:
          logging.info("RSID %s NOTFOUND (2)", rsid)
          pass
      else:
        logging.info("Concordancy check fail for - %s" % (rsid))

  #chi_test_count, allele_disc_count, overlap_count, cr_check_count = mm.get_counts()
  #logging.info("Overlap check count = %d, cr_check_count = %d", overlap_count, cr_check_count)
  chi_test_count, allele_disc_count, overlap_count = mm.get_counts()
  logging.info("CHI test count = %d, Allele discord count = %d, Overlap check count = %d", 
    chi_test_count, allele_disc_count, overlap_count)

  return count 


# execution flow starts here
#
parser = OptionParser()
parser.add_option("-r", "--rsids", dest="rsids",
   help="dbSNP rsid", metavar="STR")
parser.add_option("-f", "--snpfile", dest="snpfile",
   help="file of snps", metavar="STR")
parser.add_option("-c", "--check", dest="check",
   help="check merge candidates for concordance", metavar="STR")
parser.add_option("-p", "--prfx", dest="prfx",
   help="VCF file data path prefix", metavar="STR")
parser.add_option("-s", "--chipval", dest="chipval",
   help="min allowed chi sq p val from genotype count check", metavar="FLOAT")
parser.add_option("-l", "--logfile", dest="logfile",
   help="logfile", metavar="STR")

(options, args) = parser.parse_args()
logging.basicConfig(filename=options.logfile,level=logging.DEBUG)

if options.check == None:
  options.check = "N"

if options.chipval == None:
  options.chipval = 0.00000001
else:
  options.chipval = float(options.chipval)

if (options.rsids == None) and (options.snpfile == None):
  logging.info("Must supply either --rsids or --snpfile")
  sys.exit(0)

rec_count = main(options)
logging.info("END: %f seconds output=%d" % (time.time() - start_time, rec_count))

