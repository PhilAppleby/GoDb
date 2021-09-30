# ----------------------------------------------------------------
# Command line args:
# file1, file2, logfile
# ----------------------------------------------------------------
# 
import time
import logging
import re
import os, sys
from optparse import OptionParser
from vcfrecord import VCFrecord

start_time = time.time()

# 1st 3 are dicts of ints vs strings
file1_positions = {}
file2_positions = {}
combined_positions = {}
combined_samples = {} 
# for max comparison
max_key = 999999999999999

def load_sample_positions(fh1, fh2, vcff):
  cmts1, f1_hdr_rec = get_comments(fh1)
  cmts2, f2_hdr_rec = get_comments(fh2)

  for cmt in cmts1:
    print cmt

  f1_hdr_data = vcff.get_data_array(f1_hdr_rec)
  f2_hdr_data = vcff.get_data_array(f2_hdr_rec)
  load_f1_samples(f1_hdr_data, vcff)
  load_f2_samples(f2_hdr_data, vcff)
  return f1_hdr_data, f2_hdr_data

def get_comments(fh):
  get_cmts = True
  cmts = []
  while get_cmts == True:
    line = fh.readline().strip() 
    if line.startswith("##"):
      cmts.append(line)
    elif line.startswith("#"):
      hdr = line
      get_cmts = False
    else:
# shouldn't happen
      get_cmts = False

  logging.info("Num headers: %d" % (len(cmts)))
  return (cmts, hdr) 

def load_f1_samples(hdr_data, vcff):
  prefx, suffx = vcff.get_prfx_sfx_from_array(hdr_data) 
  for i, sample in enumerate(suffx):
    file1_positions[i] = sample
    if sample not in combined_samples:
      combined_samples[sample] = 1

def load_f2_samples(hdr_data, vcff):
  prefx, suffx = vcff.get_prfx_sfx_from_array(hdr_data) 
  for i, sample in enumerate(suffx):
    file2_positions[i] = sample
    if sample not in combined_samples:
      combined_samples[sample] = 1

def output_combined_hdr(data1, data2, vcff):
  rec_prefx = "" 
  rec_suffx = ["."] * len(combined_samples)
  rec_prefx, line_suffx = vcff.get_prfx_sfx_from_array(data1)
  for i, prochi in enumerate(line_suffx):
    sid = file1_positions[i]
    oi = combined_positions[sid]
    rec_suffx[oi] = prochi
  rec_prefx, line_suffx = vcff.get_prfx_sfx_from_array(data2)
  for i, prochi in enumerate(line_suffx):
    sid = file2_positions[i]
    oi = combined_positions[sid]
    rec_suffx[oi] = prochi
  if rec_prefx != "":
    print '\t'.join(rec_prefx + rec_suffx)

def output_combined_record(data1, data2, vcff, idxs, callrate1=0.0, callrate2=0.0):
  rec_prefx = "" 
  rec_suffx = ["."] * len(combined_samples)
  if len(data1) != 0:
    rec_prefx, line_suffx = vcff.get_prfx_sfx_from_array(data1)
    for i, genotype in enumerate(line_suffx):
      sid = file1_positions[i]
      oi = combined_positions[sid]
      rec_suffx[oi] = call_genotype(genotype, rec_suffx[oi], i, rec_prefx[2], idxs, callrate1, callrate2)
  if len(data2) != 0:
    rec_prefx, line_suffx = vcff.get_prfx_sfx_from_array(data2)
    for i, genotype in enumerate(line_suffx):
      sid = file2_positions[i]
      oi = combined_positions[sid]
      rec_suffx[oi] = call_genotype(genotype, rec_suffx[oi], i, rec_prefx[2], idxs, callrate1, callrate2)
  if rec_prefx != "":
    print '\t'.join(rec_prefx + rec_suffx)

def call_genotype(g1, g2, idx_posn, rsid, idxs, callrate1, callrate2):
  if g2 == ".":
    return g1
  if g1 == ".":
    return g2
  max1 = get_max_prob(g1, idxs)
  max2 = get_max_prob(g2, idxs)
  if max1 > max2:
    return g1
  elif max2 > max1:
    return g2
  # on equality we need to check callrates
  else: 
    if g1 != g2:
      if callrate2 > callrate1:
        logging.info("Gen mismatch CR2 GT CR1 %s, %d, %s, %.5f, %s, %.5f" % (rsid, idx_posn, g1, callrate1, g2, callrate2))
        return g2
      else:
        # favours g1 for equal callrates
        logging.info("Gen mismatch CR1 GTE CR2 %s, %d, %s, %.5f, %s, %.5f" % (rsid, idx_posn, g1, callrate1, g2, callrate2))
        return g1
  # reached when genotypes are equal
  return g1

def get_max_prob(geno, gidxs):
  #print geno
  gen_vals = geno.split(":")
  probs = gen_prob[gidxs["GP"]].split(",")
  max_prob = 0.0
  for prob in probs:
    if float(prob) > max_prob:
      max_prob = float(prob)
  return max_prob

def main(options):
  #print options.file1
  #print options.file2

  try:
    fh1 = open(options.file1, "r") 
    #fh2 = open(options.file2, "r") 
    fh2 = sys.stdin
  except IOError as e:
    logging.info("I/O error({0}): {1}".format(e.errno, e.strerror))
    exit()
  except TypeError as e:
    logging.info("Missing arguments " + e)
    exit()
  except:
    logging.info("Unexpected error:" + sys.exc_info()[0])
    exit()

  vcff = VCFrecord()

  f1_hdr, f2_hdr = load_sample_positions(fh1, fh2, vcff)
  #print len(file1_positions)
  #print len(file2_positions)

  srtd_samples = sorted(combined_samples)
  #print len(srtd_samples)
  for i, sample in enumerate(srtd_samples):
    combined_positions[sample] = i

  output_combined_hdr(f1_hdr, f2_hdr, vcff)
  line1 = fh1.readline().strip()
  data1 = vcff.get_data_array(line1)
  key1 = vcff.get_posn_from_array_as_int(data1)
  fmts1 = vcff.get_fmts_from_array(data1)
  idxs1 = vcff.get_fmt_indices(fmts1, ["GT","GP"])
  line2 = fh2.readline().strip()
  data2 = vcff.get_data_array(line2)
  key2 = vcff.get_posn_from_array_as_int(data2)
  fmts2 = vcff.get_fmts_from_array(data1)
  idxs2 = vcff.get_fmt_indices(fmts2, ["GT","GP"])
  #print key1, key2
  f1_count = 1
  f2_count = 1
  out_count = 0
  discord_count = 0

  while True:
    if (key1 == max_key and key2 == max_key):
      break
    if (key1 > key2):
      output_combined_record([], data2, vcffi, idxs1)
      out_count += 1
      line2 = fh2.readline().strip()
      if line2 == "":
        key2 = max_key
      else:
        f2_count += 1
        data2 = vcff.get_data_array(line2)
        key2 = vcff.get_posn_from_array_as_int(data2)
    elif (key2 > key1):
      output_combined_record(data1, [], vcff, idxs1)
      out_count += 1
      line1 = fh1.readline().strip()
      if line1 == "":
        key1 = max_key
      else:
        f1_count += 1
        data1 = vcff.get_data_array(line1)
        key1 = vcff.get_posn_from_array_as_int(data1)
    else:
      # On equality - check for allele concordance
      AlleleA1, AlleleB1 = vcff.get_alleles_from_array(data1)
      AlleleA2, AlleleB2 = vcff.get_alleles_from_array(data2)
      # TODO HWE concordance check - but how do we set thresholds?
      if (AlleleA1 == AlleleA2) and (AlleleB1 == AlleleB2):
        output_combined_record(data1, data2, vcff, idxs1, vcff.get_call_rate_from_array(data1), vcff.get_call_rate_from_array(data2))
        out_count += 1
      else:
        discord_count += 1
      line1 = fh1.readline().strip()
      if line1 == "":
        key1 = max_key
      else:
        f1_count += 1
        data1 = vcff.get_data_array(line1)
        key1 = vcff.get_posn_from_array_as_int(data1)
      line2 = fh2.readline().strip()
      if line2 == "":
        key2 = max_key
      else:
        f2_count += 1
        data2 = vcff.get_data_array(line2)
        key2 = vcff.get_posn_from_array_as_int(data2)
        
  fh1.close()
  fh2.close()
  return f1_count, f2_count, out_count, discord_count 


# execution flow starts here
#
parser = OptionParser()
parser.add_option("-f", "--file1", dest="file1",
   help="merge file1", metavar="STR")
parser.add_option("-a", "--ap1", dest="ap1",
   help="assay platform 1", metavar="STR")
parser.add_option("-p", "--ap2", dest="ap2",
   help="assay platform 2", metavar="STR")
#parser.add_option("-g", "--file2", dest="file2",
#   help="merge file2", metavar="STR")
parser.add_option("-l", "--logfile", dest="logfile",
   help="logfile", metavar="STR")

(options, args) = parser.parse_args()
logging.basicConfig(filename=options.logfile,level=logging.DEBUG)

f1_count, f2_count, out_count, discord_count = main(options)
logging.info("END: %f seconds f1read=%d, f2read=%d, output=%d, discord_count=%d" % (time.time() - start_time, f1_count, f2_count, out_count, discord_count))

