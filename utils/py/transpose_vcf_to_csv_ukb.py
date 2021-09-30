# Transpose a VCF file - can't handle more than a few records at a time (currently once we get above 20
# we call a halt and just produce results).
import time
import logging
import re
import os, sys
from optparse import OptionParser
from vcfrecord import VCFrecord
from mafhelper import Mafhelper # minor allele freq helper method

start_time = time.time()
assay_expand = {}
assay_expand["A"] = "affy"
assay_expand["I"] = "illumina"
assay_expand["E"] = "exome"
assay_expand["M"] = "metabo"
assay_expand["A1"] = "affy1KG"
assay_expand["I1"] = "illumina1KG"
assay_expand["B"] = "broad"
scalls = ["0/0", "0/1", "1/1", "./."]
icalls = {
    "0/0": 0,
    "0/1": 1,
    "1/0": 1,
    "1/1": 2,
    "./.": -9,
    "": ""
    }

def get_max_prob(gen_vals, probidx):
  probs = gen_vals[probidx].split(",")
  max_prob = 0.0
  max_idx = 3

  for idx, prob in enumerate(probs):
    if float(prob) > max_prob:
      max_prob = float(prob)
      max_idx = idx

  return max_prob, max_idx

def main(options):
  hdrData = ["id"]
  sampleDict = {}
  colPosns = {}
  RefAlleleDict = {}
  AltAlleleDict = {}
  count = 0

  mafh = Mafhelper()

  for line in sys.stdin:
    line = line.strip()
    if (line.startswith('##')):
      pass
    else:
      vcfr = VCFrecord(line)
      prfx, sfx = vcfr.get_prfx_sfx()
      #print prfx
      if (line.startswith('#')):
# Parse out the header record.
        for i, col_hdr in enumerate(sfx):
          colPosns[i] = col_hdr
          sampleDict[col_hdr] = []
      else:
        flip = False
        varid = vcfr.get_varid_ukb()
        #logging.info("varid=%s", varid)
        ref, alt = vcfr.get_alleles()
        probidx = vcfr.get_probidx()
        hdr_allele = alt
        homref_count, het_count, homalt_count, nc_count, miss_count = vcfr.get_allele_counts()
        call_count = homref_count + het_count + homalt_count
        maf, ma = mafh.maf(het_count, homref_count, ref, homalt_count, alt, nc_count)
        RefAlleleDict[varid] = ref
        AltAlleleDict[varid] = alt
        #if ma == ref:
        #  flip = True
        #  hdr_allele = ref
        #  logging.info("FLIP for %s, %s, %s", varid, ref, alt)
        hdrData.append(varid)
        for i, str_geno in enumerate(sfx):
          if str_geno != ".":
            geno = str_geno.split(":")
            max_prob, max_idx = get_max_prob(geno, probidx)
            i_call = icalls[geno[0]]
            if flip == True:
              if i_call == "0":
                i_call == "2"
              elif i_call == "2":
                i_call = "0"
            sampleDict[colPosns[i]].append(str(i_call))
          else:
            sampleDict[colPosns[i]].append("")

  print ",".join(hdrData)
  for samp in sampleDict:
    count += 1
    print ",".join([samp] + sampleDict[samp])
  return count


# execution flow starts here
#
parser = OptionParser()
parser.add_option("-l", "--logfile", dest="logfile",
    help="logfile", metavar="STR")
parser.add_option("-p", "--platform", dest="platform",
    help="Optional assay platform", metavar="STR")

(options, args) = parser.parse_args()
logging.basicConfig(filename=options.logfile,level=logging.DEBUG)

rec_count = main(options)
logging.info("END: %f seconds output=%d" % (time.time() - start_time, rec_count))
