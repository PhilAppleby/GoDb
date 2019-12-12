import pymongo
import time
import logging
from zipfile import ZipFile
from zipfile import ZIP_DEFLATED
from optparse import OptionParser

from config import probidx

from models_test import DataStore


# execution flow starts here
#
parser = OptionParser()
parser.add_option("-f", "--filename", dest="filename",
    help="rsid file", metavar="STR")
(options, args) = parser.parse_args()

connection = pymongo.MongoClient("mongodb://localhost")
db = connection.genomicsdb_test
#db = connection.gwas_metadata

ds = DataStore(db, "genomicsdb_test", "akh", False, probidx)
start_time = time.time()
calls = ["0/0", "0/1", "1/1", "Missing"]
icalls = [0, 1, 2, -9]


def get_call(probs, threshold):
  max_prob = 0.0 
  max_idx = 3 

  #print "Probs:", probs
  for idx, prob in enumerate(probs):
    if float(prob) > max_prob:
      max_prob = float(prob)
      max_idx = idx 

  if (threshold !=0.0):
    #print 'threshold', threshold, max_prob
    if max_prob < threshold:
      #print 'LT threshold', threshold, max_prob
      max_idx = 3

  return (calls[max_idx], icalls[max_idx], max_prob)


#print ds.get_prochi_from_mprochi("CAO4147247")
def resolve_geno(genlist, threshold):
  maxprob = 0.0 
  maxidx = -1
  if len(genlist) == 1:
    return genlist[0]
  elif len(genlist) > 1:
    for idx, gendata in enumerate(genlist):
      if gendata[1] == 0:
        #print "Decided on D Type:", genlist[idx][0], rsid, samp
        return genlist[idx]
      dataVals = gendata[2].split(':')
      probVals = dataVals[1].split(',')
      (probcall, intcall, outprob) = get_call(probVals, threshold)
      if outprob > maxprob:
        maxprob = outprob
        maxidx = idx 
  if maxprob > 0.0:
    #print "Decided on prob:", maxprob, maxidx, rsid, samp
    return genlist[maxidx]
  return []
  
#print ds.get_prochi_from_plateid("75246_D09")

#samples = ds.get_converted_samples()
#for sample in samples:
#  print sample
#totals = ds.get_marker_totals()

#for vals in totals:
#  print vals[0], vals[1]

#print ds.get_sample_count()
#print ds.get_all_samples()
#(lres, snpres, msg) = ds.get_rslist_data("/home/hicadmin/data/txt", "moneeza_2.txt", 0.9, "csv_calls")
#(lres, snpres, msg) = ds.get_rslist_data("/home/hicadmin/data/txt", "pda_rsid.txt", 0.9, "csv_calls")
print "START"
logging.basicConfig(filename="./tmp.log",level=logging.DEBUG)
count = 0
(pres, snpres, msg) = ds.get_rslist_file_data(options.filename, 0.9, [])
#(pres, snpres, msg) = ds.get_rslist_data(["rs7412"], 0.9, [])

ares = {}

for assaytype in pres:
  ares[assaytype] = '\n'.join(pres[assaytype])
  
connection.close()

print "END:", time.time() - start_time, "seconds", count
zipname = "/homes/hic-gendb/prod/data/results/result.zip"

with ZipFile(zipname, 'w') as resZip:
  resZip.writestr('snp_data.csv', snpres, ZIP_DEFLATED)
  #resZip.writestr('combined_sample_data.csv', res, ZIP_DEFLATED)
  for assaytype in ares:
    resZip.writestr(assaytype + '_sample_data.csv', ares[assaytype], ZIP_DEFLATED)


