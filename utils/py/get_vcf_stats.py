# ----------------------------------------------------------------
# Command line args:
# ----------------------------------------------------------------
#
import time
import logging
import re
import os, sys
from optparse import OptionParser

from vcfrecord import VCFrecord # Record access via VCF helper methods
from mafhelper import Mafhelper # minor allele freq helper method
from hwehelper import Hwehelper # Hardy Weinberg Equilibrium Fischer Exact Test helper method

start_time = time.time()

def main():
  mafh = Mafhelper()
  hweh = Hwehelper()
  in_count = 0
  hdr_count = 0
  homr_total = 0
  het_total = 0
  homa_total = 0
  virt_nc_total = 0
  miss_total = 0

  print "SNPId,AssayType,chr,pos,REF,ALT,Minor,MAF,CallRate,HWE_pval"

  for line in sys.stdin:
    line = line.strip()
    in_count += 1
    if line.startswith("#"):
      hdr_count += 1
      continue


    vcfr = VCFrecord(line)
    varid = vcfr.get_varid_ukb()
    chromosome = vcfr.get_chr()
    posn = vcfr.get_posn_as_int()
    ref, alt = vcfr.get_alleles()
    homref_count, het_count, homalt_count, virt_nc_count, miss_count = vcfr.get_allele_counts()
    call_count = homref_count + het_count + homalt_count
    #nocall_count = virt_nc_count + miss_count
    nocall_count = virt_nc_count
    call_rate = float(call_count) / float(call_count + nocall_count)
    homr_total += homref_count
    het_total += het_count
    homa_total += homalt_count
    virt_nc_total += virt_nc_count
    miss_total += miss_count
    try:
      hwe = hweh.HWE_exact(het_count, homref_count, homalt_count, call_count)
      maf, ma = mafh.maf(het_count, homref_count, ref, homalt_count, alt, virt_nc_count)
    except ZeroDivisionError:
      logging.info("DIV 0 error at %d (%d), where hom_r=%d, het=%d, home_a=%d, cc=%d", in_count, posn, homref_count, het_count, homalt_count, call_count)
    print "%s,combo,%s,%d,%s,%s,%s,%s,%.3f,%s" % (varid, chromosome, posn, ref, alt, ma, maf, call_rate, hwe)
  return in_count, hdr_count, homr_total, het_total, homa_total, virt_nc_total, miss_total


# execution flow starts here
#
#parser = OptionParser()

#(options, args) = parser.parse_args()

FORMAT = '%(asctime)-15s %(message)s'
logging.basicConfig(format=FORMAT, filename="./log.txt",level=logging.DEBUG)

in_count, hdr_count, homr_total, het_total, homa_total, virt_nc_total, miss_total = main()
logging.info("END: %f seconds rec_read=%d, hdrs=%d, homref=%d, het=%d, homalt=%d, virt_nocall=%d, miss=%d" % (time.time() - start_time, in_count, hdr_count, homr_total, het_total, homa_total, virt_nc_total, miss_total))
