#
#
import time
import datetime
import os, sys
from grshelper import GRShelper # GRS utilities

from optparse import OptionParser

# SNP statistics data is derived from local genotype data
def load_snp_local_stats_data(fh):
  alleles = {} # always return REF/ALT
  minor_allele = {} # Usually ALT but can be REF
  maf = {} # Refers to the minor allele

  hdr = fh.readline().strip().split(',')

  for line in fh:
    data=line.strip().split(',')
    alleles[data[0]] = data[4] + "/" + data[5]
    minor_allele[data[0]] = data[6]
    maf[data[0]] = float(data[7])

  return alleles, minor_allele, maf

def main(options):
  count = 0
  rejected_snps = {}
  flip_snps = {}
  # Housekeeping: open input files and load reference data into memory
  try:
    fh = open(options.snpfile, "r")
    fh2 = open(options.statsfile, "r")
    grsh = GRShelper()
    rsids, snp_wts, wt_alleles, wt_rafs = grsh.load_snp_weight_data_common(fh)
    local_alleles, local_minor_alleles, local_mafs = load_snp_local_stats_data(fh2)
  except IOError as e:
    print "I/O error({0}): {1}".format(e.errno, e.strerror)
    exit()
  except TypeError as e:
    print "Missing arguments ", e
    exit()
  except:
    print "Unexpected error: %s", str(sys.exc_info())
    sys.exit()

  print "STATS num snp_wts = %d" % (len(snp_wts))
  numsnps = len(rsids)
  cthreshold = numsnps / 2

  # Process each rsid in the snp weight data in turn
  for rsid in rsids:
    if rsid not in local_alleles:
      print "NOTFOUND %s" % (rsid)
      continue
    wtRa = wt_alleles[rsid]
    locall = local_alleles[rsid].split('/')
    locNra = locall[0]
    locRa = locall[1]
    maf = float(local_mafs[rsid])
    locma = local_minor_alleles[rsid]
    locRaf = maf
    if locma == locNra:
      locRaf = (1.0 - float(maf))

    if not grsh.isCompatible(wtRa, wt_rafs[rsid], locRa, locNra, maf):
      rejected_snps[rsid] = 1
      print "REJ: %s" % (rsid)
      continue

    #if grsh.notEqual(wtRa, locRa):
    if grsh.should_flip_full(rsid, wtRa, float(wt_rafs[rsid]), locNra, locRa, locRaf):
      print "FLIPPED for %s, wgt=%f their side: ra=%s, raf=%f, local: ra=%s, nra=%s, [ma=%s], raf=%f" % (rsid, float(snp_wts[rsid]), wtRa, float(wt_rafs[rsid]), locRa, locNra, locma, locRaf)
      flip_snps[rsid] = 1

  for snpid in snp_wts:
    locMult = 1
    if snpid in flip_snps and options.flipgeno == False:
      locMult = -1
    print "SNPWGT,%s,%f" % (snpid, float(snp_wts[snpid]) * locMult)

  snp_arr = []
  csvhdr = sys.stdin.readline().strip().split(',')
  snp_arr = csvhdr[1:]
  print "OUTPUT,id,count,score"

  for line in sys.stdin:
    data = line.strip().split(',')
    gr_score = 0.0
    snp_count = 0
    sampleid = data[0]
    for i, geno in enumerate(data[1:]):
      #print "GENOVAL,SAMPLEID,RSID,%s,%s,%s" % (sampleid, snp_arr[i], geno)
      rsid = snp_arr[i]
      mult = 1.0
      if rsid in rejected_snps:
        continue
#     Only count genotypes of value 0,1,2
      if geno == '-9':
        continue
      if geno == 'NA':
        continue
      if geno == '':
        continue
      if rsid in flip_snps:
        if options.flipgeno == True:
          #print "FLIPGENO,SAMPLEID,RSID,%s,%s,%s" % (sampleid, snp_arr[i], geno)
          if geno == '2':
            geno = '0'
          elif geno == '0':
            geno = '2'
        else:
          #print "FLIPSIGN,SAMPLEID,RSID,%s,%s,%s" % (sampleid, snp_arr[i], geno)
          mult = -1.0
      snp_count += 1
      wt = float(snp_wts[snp_arr[i]])
      wt = wt * mult
      wtra = wt_alleles[rsid]
      #print "SCORE %d * %f, %f" % (int(geno), float(snp_wts[snp_arr[i]]), (int(geno) * float(snp_wts[snp_arr[i]])))
      gr_score += int(geno)*wt
      #print "GENO %s %s %s %s*%f,[%f]" % (rsid, wtra, sampleid, str(geno), wt, int(geno)*wt)
    if snp_count >= cthreshold:
      print "OUTPUT,%s,%d,%f" % (sampleid, snp_count, gr_score)
    else:
      print "LOWCOUNT,%s,%d,%f" % (sampleid, snp_count, gr_score)

  return count


# execution flow starts here
#
start_time = time.time()
parser = OptionParser()
#
parser.add_option("-p", "--snpfile", dest="snpfile",
   help="file of snps", metavar="STR")
parser.add_option("-s", "--statsfile", dest="statsfile",
   help="file of local snp statistics", metavar="STR")
parser.add_option("-f", "--flipgeno", dest="flipgeno",
   help="Should we apply genotype flipping?", metavar="STR")

start_time = time.time()
(options, args) = parser.parse_args()

if (options.flipgeno == None) or (options.flipgeno == "N"):
  options.flipgeno = False
else:
  options.flipgeno = True

count = main(options)
print("END: {0:.5f} seconds output={1}".format(time.time() - start_time, count))
#print "END:", time.time() - start_time, "seconds", count
