# Load the mongodb metadata markers collection
#
import time
import os, sys
import csv
import logging
from optparse import OptionParser
from dbsnpfile import Dbsnpfile # File access via tabix helper methods
from vcfrecord import VCFrecord
from godb import GoDb

start_time = time.time()

def load_snp_summary(fh):
  snp_data = {}
  fh.readline()
  for line in fh:
    data = line.strip().split(',')
    snp_data[data[1]] = data[2].split(':')
  return snp_data

def get_dbsnp_rsid(dbsnpfile, chrom, posn):
  dbsnprec = dbsnpfile.get_dbsnp_file_record(options.dbsnpfile, chrom, int(posn))
  rsid = ""
  refallele = ""
  if dbsnprec != None:
    dbvcf = VCFrecord(dbsnprec)
    rsid = dbvcf.get_varid()
    refallele, altallele = dbvcf.get_alleles()
  return rsid, refallele

def get_posn_allele(snpdata):
  posndata = snpdata[1].split("_")
  allele = ""
  if len(posndata) == 2:
    allele = posndata[1]
  else:
    allele = snpdata[2]

  return posndata[0], allele

def flip_geno(geno):
  if geno == 0:
    return "2"
  if geno == 2:
    return "0"
  return str(geno)

def main(options):

  try:
    csvfile =  open(options.csvfile, "r")
    csvreader = csv.reader(csvfile)
    fh = open(options.snpsummary, "r")
    snp_data = load_snp_summary(fh)
    dbsnpfile = Dbsnpfile()
    dbsnpfile.set_tabix_file(options.dbsnpfile)
    godb = GoDb()
  except IOError as e:
    logging.fatal("I/O error({0}): {1}".format(e.errno, e.strerror))
    sys.exit()
  except TypeError as e:
    logging.fatal("Missing arguments ", e)
    sys.exit()
  except:
    logging.fatal("Unexpected error:", sys.exc_info())
    sys.exit()

  hdr = []
  hdrlen = 0
  count = 0
  flipidx = {}

  hdr = csvreader.next()
  hdrlen = len(hdr)

  for i, varid in enumerate(hdr):
    if "ID" in varid:
      pass
    elif varid.startswith("rs"):
      # Split to get the allele component, check allele va alleleA
      # If necessary add an entry to the flip array, replace the hdr element with raw rsNumber
      var = varid.split("_")
      vardata = godb.get_one_variant(var[0])
      #print var[0], var[1], vardata["alleleA"], vardata["chromosome"], vardata["position"], i
      hdr[i] = var[0]
      if var[1] == vardata["alleleA"]:
        flipidx[i] = True
    else:
      coldata = varid.split(":")
      posn, allele = get_posn_allele(coldata)
      var, ref = get_dbsnp_rsid(dbsnpfile, coldata[0], posn)
      if allele == ref:
        flipidx[i] = True
      #print coldata, var, ref, i
      hdr[i] = var

  print ",".join(hdr)

  for row in csvreader:
    count += 1
    for i, genotype in enumerate(row):
      if i in flipidx:
        row[i] = flip_geno(int(genotype))
    print ",".join(row)
  return count, hdrlen
#
# execution flow starts here
#
parser = OptionParser()
# csv file containing detail columns
parser.add_option("-c", "--csvfile", dest="csvfile",
  help="csv file containing genotype detail data", metavar="FILE")
# csv file containing all snp summaries
parser.add_option("-s", "--snpsummary", dest="snpsummary",
  help="csv file containing genotype detail data", metavar="FILE")
# tabix indexed file with dbsnp (rsid) information
parser.add_option("-d", "--dbsnpfile", dest="dbsnpfile",
   help="Path to file containing dbsnp rsids", metavar="STR")
# logging
parser.add_option("-l", "--logfile", dest="logfile",
  help="logfile", metavar="STR")

(options, args) = parser.parse_args()
logging.basicConfig(filename=options.logfile,level=logging.DEBUG)

rec_count, hdrlen = main(options)
logging.info("END: %f seconds output=%d, hdrlen=%d" % (time.time() - start_time, rec_count, hdrlen))
