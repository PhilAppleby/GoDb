import pymongo
import time
from zipfile import ZipFile
from zipfile import ZIP_DEFLATED
from optparse import OptionParser

from models import DataStore


# execution flow starts here
#
parser = OptionParser()
parser.add_option("-f", "--filename", dest="filename",
    help="rsid file", metavar="STR")
(options, args) = parser.parse_args()

connection = pymongo.MongoClient("mongodb://localhost")
db = connection.genomicsdb_test
#db = connection.gwas_metadata

ds = DataStore(db, "genomicsdb_test", "bpf", False)
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
(rslist, assaytypes, sampleDict, lres, pres, snpres, msg) = ds.get_rslist_file_data_test(options.filename, 0.9, "csv")
idx = 0
rsidx= {}
for rsid in rslist:
  rsidx[rsid] = idx
  #print rsid, idx
  idx +=1

print "SAMPLEDICT"
print "*****"
for samp in sampleDict:
  #print samp, sampleDict[samp]
  output_lines = {}
  empty_output_lines = {}
  output_lines['combined'] = ["" for x in range(len(rslist) * 4)]
  for assaytype in assaytypes:
    output_lines[assaytype] = ["" for x in range(len(rslist) * 4)]
    empty_output_lines[assaytype] = True
  for rsid in sampleDict[samp]:
    if len(sampleDict[samp][rsid]) > 0: # if a sample wasn't genotyped on a particular platform there might not be data
      idxoffset = rsidx[rsid] * 4
      geno_data = resolve_geno(sampleDict[samp][rsid], 0.9)
      dataVals = geno_data[2].split(':')
      probVals = dataVals[1].split(',')
      (probcall, intcall, outprob) = get_call(probVals, 0.9)
      output_lines['combined'][idxoffset] = str(intcall)
      output_lines['combined'][idxoffset+1] = str(outprob)
      output_lines['combined'][idxoffset+2] = geno_data[0]
      output_lines['combined'][idxoffset+3] = geno_data[4]
      for geno_data in sampleDict[samp][rsid]:
        dataVals = geno_data[2].split(':')
        probVals = dataVals[1].split(',')
        (probcall, intcall, outprob) = get_call(probVals, 0.9)
        output_lines[geno_data[0]][idxoffset] = str(intcall)
        output_lines[geno_data[0]][idxoffset+1] = str(outprob)
        output_lines[geno_data[0]][idxoffset+2] = geno_data[0]
        output_lines[geno_data[0]][idxoffset+3] = geno_data[4]
        empty_output_lines[geno_data[0]] = False

  print samp, output_lines
  print samp, empty_output_lines
#print lres
count = 0
#for r in lres:
#  res += r
#  count += 1
#  print r
connection.close()

print "END:", time.time() - start_time, "seconds", count
#zipname = "/home/hicadmin/data/results/result.zip"

#with ZipFile(zipname, 'w') as resZip:
#  resZip.writestr('sample_data.csv', res, ZIP_DEFLATED)
#  resZip.writestr('snp_data.csv', snpres, ZIP_DEFLATED)
