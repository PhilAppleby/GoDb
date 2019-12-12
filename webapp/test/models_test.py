import pymongo
import logging
import sys
import pysam
from hwe import HWE_exact
from maf import maf
from maf import af
from gwasdb import Gwasdb
from vcfrecord import VCFrecord
from zipfile import ZipFile
from zipfile import ZIP_DEFLATED

class DataStore():
  def __init__(self, db, dbname, projpref="akh", get_anochi=False, probidx=1):
    #print "Anochi logical", get_anochi
    self.db = db
    self.dbname = dbname
    self.gwasdb = Gwasdb(db)
    self.filedata_coll = _filedata(db)
    self.sam_coll = _samples(db)
    self.mkr_coll = _markers(db, self.filedata_coll, self.sam_coll, self.gwasdb, probidx)
    self.prochi_coll = _prochi_map(db, projpref, get_anochi)
    self.marker_totals = []
    self.sample_count = -1
    self.call_rates = {}
    self.probidx = probidx
    self.vcfr = VCFrecord()

  def get_db_name(self):
    return self.dbname

  def get_probidx(self):
    return self.probidx

  def make_selection_key(self, varid, assaytype):
    return varid + "_" + assaytype

  def get_rsid_prochi_data(self, rsid, prochi, threshold, filefmt):
    """
      Get  data for the prochi, dict {platform:position_in_list}
      Get marker data for the rsid (up to num of platforms)
    """
    assaytype_list_posns = {}

    sdocs = self.sam_coll.get_sampledata(prochi)

    for sdoc in sdocs:
      #print sdoc["assaytype"]
      assaytype_list_posns[sdoc["assaytype"]] = sdoc["list_posn"]

    genotypes = self.mkr_coll.get_geno_data(rsid, prochi, assaytype_list_posns)

    for genotype in genotypes:
      print genotype, genotypes[genotype]

  def get_marker_data_for_file(self, filepath, threshold):
    msg = None
    try:
      f = open(filepath, "r")
    except IOError as e:
      msg = filepath + ":" + e.strerror
      return([],[],msg)
    count = 0
    
    rslist = []
    for line in f:
      count += 1
      if count > 500:
        msg = "line count for %s gt the limit (%d)" % (filepath, 500)  
        return ([], [], msg)
      line = line.strip()
      elems = line.split()
      rslist.append(elems[0])

    f.close()
    marker_data = []
    msg = ""
    for rsid in rslist:
      (marker_docs, tmsg) = self.get_marker_summary_probs(rsid, threshold)
      if (tmsg != ""):
        msg+=tmsg
      for doc in marker_docs:
        marker_data.append(doc)
    return marker_data, msg

  def get_marker_data_by_range(self, chr, start, end, threshold=0.9):
    rslist = []
    marker_data = []
    msg = ""
    return(self.mkr_coll.get_marker_data_by_range(chr, start, end))

  def get_range_data(self, chr, start, end, threshold, download_list):
    (docs, msg) = self.mkr_coll.get_marker_data_by_range(chr, start, end)
    if len(docs) == 0:
      return([], [], msg)

    rsdict={}
    rslist=[]
    for doc in docs:
      #print"RANGE rsid", doc["rsid"], doc["assaytype"], doc["position"]
      rsdict[doc["rsid"]] = 1

    for rsid in rsdict:
      rslist.append(rsid)

    #return([], [], "")
    #print "Call get rslist data"
    return(self.get_rslist_data(rslist, threshold, download_list))

  def build_csv_data(self, rslist, sampleDict, assaytypes, threshold):
    # NOTE: maintaining rslist order is vital!
    rsidx= {}
    normalised = False
    idx = 0
    for rsid in rslist:
      rsidx[rsid] = idx
      idx +=1
    
    assaytypes['combined'] = 1
    by_platform_data = {}
    for assaytype in assaytypes:
      by_platform_data[assaytype] = []
    hdrData = ["sampleId"]
    for rsid in rslist:
      #print "RSID", rsid
      hdrData.append(rsid)
      hdrData.append(rsid + "_c")
      hdrData.append(rsid + "_p")
      hdrData.append(rsid + "_alt")
    hdrString = ','.join(hdrData)
    
    for assaytype in assaytypes:
      by_platform_data[assaytype].append(hdrString)
    for samp in sampleDict:
      output_lines = {}
      filled_output_lines = {}
      filled_output_lines['combined'] = True
      for assaytype in assaytypes:
        output_lines[assaytype] = ["" for x in range(len(rslist) * 4)] 
      for rsid in sampleDict[samp]:
        if len(sampleDict[samp][rsid]) > 0: # if a sample wasn't genotyped on any platform there might not be data
          idxoffset = rsidx[rsid] * 4 
          # resolve_geno is at the crux - need to change to test CR?
          #logging.info("Call resolve_geno %s, %s", samp, str(sampleDict[samp][rsid]))
          geno_data = self.resolve_geno(sampleDict[samp][rsid], rsid, samp, threshold)
          dataVals = geno_data[2].split(':')
          probVals = dataVals[self.get_probidx()].split(',')
          (probcall, intcall, outprob) = self.mkr_coll.get_call(probVals, threshold)
          #intcall, normalised = self.get_integer_call(intcall, geno_data[5], geno_data[6])
          output_lines['combined'][idxoffset] = str(intcall)
          output_lines['combined'][idxoffset+1] = str(outprob)
          output_lines['combined'][idxoffset+2] = geno_data[0]
          output_lines['combined'][idxoffset+3] = geno_data[4]
          for geno_data in sampleDict[samp][rsid]:
            dataVals = geno_data[2].split(':')
            probVals = dataVals[self.get_probidx()].split(',')
            (probcall, intcall, outprob) = self.mkr_coll.get_call(probVals, threshold)
            #intcall, normalised = self.get_integer_call(intcall, geno_data[5], geno_data[6])
            output_lines[geno_data[0]][idxoffset] = str(intcall)
            output_lines[geno_data[0]][idxoffset+1] = str(outprob)
            output_lines[geno_data[0]][idxoffset+2] = geno_data[0]
            output_lines[geno_data[0]][idxoffset+3] = geno_data[4]
            filled_output_lines[geno_data[0]] = True
      
      for assaytype in assaytypes:
        if assaytype in filled_output_lines:
          by_platform_data[assaytype].append(samp + "," + ",".join(output_lines[assaytype]))
    return by_platform_data
  
  def build_gen_data(self, rslist, sampleDict, threshold):
    data = []
    hdrData = ["rsid"]
    for sampleid in sorted(sampleDict):
      hdrData.append(sampleid)
      hdrData.append(sampleid)
      hdrData.append(sampleid)
    hdrString = ' '.join(hdrData)
    data.append(hdrString)
    for rsid in rslist:
      line = rsid + " "
      for samp in sorted(sampleDict):
        geno_data = self.resolve_geno(sampleDict[samp][rsid], rsid, samp, threshold)
        if len(geno_data) > 2:
          dataVals = geno_data[2].split(':')
          probVals = dataVals[self.get_probidx()].split(',')
          line += ' '.join(probVals) + " "
        else:
          line += ' '.join(["0","0","0"]) + " "
        
      data.append(line[:-1])
    return data
  
  def get_integer_call(self, intcall, afreq, bfreq, normalise=False):
    rtncall = intcall
    normalised = False
    #print rtncall, normalise, afreq, bfreq, 
    if normalise == True:
      #print "normalising",
      if afreq < bfreq:
        normalised = True
        #print "LT",
        if rtncall == 2:
          #print "2",
          rtncall = 0
        elif rtncall == 0:
          #print "0",
          rtncall = 2
        #else:
          #print "1",
      #else:
        #print "GE",
    #print "return", rtncall
    return rtncall, normalised

  def resolve_geno(self, genlist, rsid, samp, threshold):
    maxprob = 0.0
    maxidx = -1
    if len(genlist) == 1:
      return genlist[0]
    elif len(genlist) > 1:
      for idx, gendata in enumerate(genlist):
        #if gendata[1] == 0:
        #  #print "Decided on D Type:", genlist[idx][0], rsid, samp
        #  return genlist[idx]
        dataVals = gendata[2].split(':')
        probVals = dataVals[self.get_probidx()].split(',')
        (probcall, intcall, outprob) = self.mkr_coll.get_call(probVals, threshold)
        if outprob > maxprob:
          maxprob = outprob
          maxidx = idx
    if maxprob > 0.0:
      #print "Decided on prob:", maxprob, maxidx, rsid, samp
      return genlist[maxidx]
    return []

  def get_marker_summary_probs(self, rsid, threshold):
    marker_array = []
    msg = ""
    #print "get_marker_summary_probs:", threshold
    docs = self.mkr_coll.get_marker_data_multi(rsid)
    for doc in docs:
      fpath = self.filedata_coll.get_filepath(doc["assaytype"], doc['chromosome'])
      # this is not ideal - need to get on top of this chromosome id thing
      chr = "%.2d" % int(doc["chromosome"])
      rec, fullrec = self.mkr_coll.get_raw_marker_values(fpath, doc["rsid"], chr, doc['position'])
      #print "REC", rec
      prfx, sfx = self.vcfr.get_prfx_sfx_from_array(rec)
      (gc_count_dict, sample_count, geno_count, maf, alleleAf, alleleBf, p_hwe) = self.mkr_coll.get_genotype_probs(sfx, threshold)
      doc['selected'] = 1
      doc['a_af'] = alleleAf
      doc['b_af'] = alleleBf
      doc['hwe_p'] = p_hwe
      doc['Missing'] = 0
      if 'Missing' in gc_count_dict:
        doc['Missing'] = gc_count_dict['Missing']
      marker_array.append(doc)
    if len(marker_array) == 0:
      msg = "Variant NOT FOUND - %s, " % (rsid)
    return (marker_array, msg)

  def get_marker_data(self, markerid, assaytype):
    return self.mkr_coll.get_marker_data(markerid, assaytype)

  def get_sample(self, sampleid):
    return self.sam_coll.get_sample(sampleid)

  def get_marker_totals(self):
    if len(self.marker_totals) == 0:
      self.marker_totals = self.mkr_coll.get_marker_totals()
    return self.marker_totals

  def get_sample_count(self, assaytype):
    if self.sample_count == -1:
      self.sample_count = self.sam_coll.get_count(assaytype)
    return self.sample_count
  
  def get_all_samples(self):
    return self.sam_coll.get_all_samples()
  
  def convert_to_prochi(self, cvt_value):
    """
    Convert the supplied value to prochi by whichever method works (or return supplied value if all else fails)
    """
    rtn_value = self.get_prochi_from_mprochi(cvt_value)
    if rtn_value == None:
      rtn_value = self.get_prochi_from_plateid(cvt_value)
    if rtn_value == None:
      rtn_value = cvt_value
    return rtn_value

  def get_prochi_from_mprochi(self, mprochi):
    """
    Get the prochi_maps value for the supplied arg
    """
    return self.prochi_coll.get_anochi_or_prochi_from_mprochi(mprochi)

  def get_prochi_from_plateid(self, plateid):
    """
    Get the prochi_maps value for the supplied arg
    """
    return self.prochi_coll.get_prochi_from_plateid(plateid)
  
  def get_converted_samples(self):
    return [self.convert_to_prochi(samp) for samp in self.sam_coll.get_all_samples()]

  def make_zipfile(self, sample_return_data, snp_return_data, uploadDir, zipfilename):
    """
    moved here from views.py
    """
    ares = {}
    for assaytype in sample_return_data:
      ares[assaytype] = '\n'.join(sample_return_data[assaytype])
    zipname = uploadDir + "/" + zipfilename
    with ZipFile(zipname, 'w') as resZip:
      resZip.writestr('snp_summary.csv', snp_return_data, ZIP_DEFLATED)
      for assaytype in ares:
        resZip.writestr(assaytype + '_samples.csv', ares[assaytype], ZIP_DEFLATED)
    with open(zipname, 'r') as f:
      body = f.read()
    return(body)
    response = make_response(body)
    response.headers["Content-Disposition"] = "attachment; filename=" + zipfilename
    return(response)
  
  def get_rslist_file_data(self, filepath, threshold, download_list):
    msg = None
    try:
      f = open(filepath, "r")
    except IOError as e:
      msg = filepath + ":" + e.strerror
      return([],[],msg)
    count = 0
    
    rslist = []
    for line in f:
      count += 1
      if count > 500:
        msg = "line count for %s gt the limit (%d)" % (filepath, 500)  
        return ([], [], msg)
      line = line.strip()
      elems = line.split()
      rslist.append(elems[0])

    f.close()
    logging.info("get_rslist_file_data: %.2f", float(threshold))
    return(self.get_rslist_data(rslist, threshold, download_list))

  def get_rslist_data(self, input_rslist, threshold, download_list):
    msg = None
    snpdata = "SNPId,AssayType,chr,pos,REF,REF_fr,ALT,ALT_fr,MAF,Imputed,CallRate,HWE_pval,Info\n"

    data = []
    rslist = []
    assaytypes = {}
    Afreq = {}
    Bfreq = {}

    data_count = 0
    impDict = {}
    for rsid in input_rslist:
      docs = self.mkr_coll.get_marker_data_multi(rsid)
      if len(docs) > 0:
        rslist.append(rsid)
      # handling SNPs on multiple platforms
      for doc in docs:
        #print doc
        # first get filepath
        #select_key = self.make_selection_key(doc["rsid"], doc["assaytype"])
        #if select_key in download_list:
        fpath = self.filedata_coll.get_filepath(doc["assaytype"], doc["chromosome"])
        # this is not ideal - need to get on top of this chromosome id thing
        chr = "%.2d" % int(doc["chromosome"])
        # get raw marker data
        rec, fullrec = self.mkr_coll.get_raw_marker_values(fpath, doc["rsid"], chr, doc['position'])
        geno_count = 0
        sample_count = 0
        hwep = 0.0
        #print "REC", rec, fpath
        prfx, sfx = self.vcfr.get_prfx_sfx_from_array(rec)
        (gc_count_dict, sample_count, geno_count, maf, alleleAf, alleleBf, p_hwe) = self.mkr_coll.get_genotype_probs(sfx, threshold)
        Afreq[doc["rsid"] + "_" + doc["assaytype"]] = alleleAf
        Bfreq[doc["rsid"] + "_" + doc["assaytype"]] = alleleBf
        data_count += 1
        hwep = float(p_hwe) 
          
        # add record to list of all records for the rslist
        data.append(doc["assaytype"] + '\t' + fullrec)
        if doc["assaytype"] not in assaytypes:
          assaytypes[doc["assaytype"]] = 1

        imputed = 0
        if "imputed" in doc:
          imputed = 1
        impDict[doc["rsid"] + "_" + doc["assaytype"]] = imputed
        #print doc["rsid"], doc["assaytype"]
        if "cohort_1_hwe" in doc:
          hwep = doc["cohort_1_hwe"]
        snpdata += "%s,%s,%s,%d,%s,%s,%s,%s,%s,%d,%.5f,%.5f,%.5f\n" % (doc["rsid"], doc["assaytype"], doc["chromosome"], doc["position"], doc["alleleA"], alleleAf, doc["alleleB"], alleleBf, maf, imputed, float(geno_count) / sample_count, hwep, doc["info"])
        #else: # not selected
          #logging.info("Line was unselected: %s", select_key)

    pdata = self.get_sample_values(data, data_count, rslist, impDict, assaytypes, Afreq, Bfreq, threshold)
    return(pdata, snpdata, msg)

  def get_sample_values(self, records, numrecs, rslist, impDict, platforms, Afreq, Bfreq, threshold):
    """Process vcf records
      NOTE: impDict keyed on a composite of rsid_platform 
       """
    #print 'START 1', len(records)
    #print impDict
    first_sample_idx = 9 + 1 # add 1 due to forcing assaytype in as col 0
    samplesByAt = {}
    for platform in platforms:
      samplesByAt[platform] = [self.convert_to_prochi(samp) for samp in self.sam_coll.get_samples(platform)]
    # A dict of dicts of tables
    sampleDict = {}
    dupDict = {}
    dupcount = 0
    for platform in samplesByAt:
      for samp in samplesByAt[platform]:
        if samp not in sampleDict:
          sampleDict[samp] = {}
          for rsid in rslist:
            sampleDict[samp][rsid] = []

    #print '2'
    count = 0
    rscount = 0
    values_totals = 0
    for line in records:
      linedata = line.split('\t')
      platform = linedata[0]
      chr = linedata[1]
      pos = linedata[2]
      rsid = linedata[3]
      alleleA = linedata[4]
      alleleB = linedata[5]
      #print "rec", rscount, rsid, platform, pos, len(linedata)
      if platform not in samplesByAt:
        print platform, "not cached"
        samplesByAt[platform] = [self.convert_to_prochi(samp) for samp in self.sam_coll.get_samples(platform)]
      for idx, elem in enumerate(linedata):
        #print idx,
        if idx >= first_sample_idx:
          arridx = idx - first_sample_idx
          if samplesByAt[platform][arridx] in sampleDict:
            sampleId = samplesByAt[platform][arridx]
            values_totals += 1
            sampleDict[sampleId][rsid].append([platform,impDict[rsid  + "_" + platform],linedata[idx],alleleA,alleleB, 
              Afreq[rsid  + "_" + platform], Bfreq[rsid  + "_" + platform]])
      print "x"
      rscount += 1 

    print '3'
    by_platform_data = self.build_csv_data(rslist, sampleDict, platforms, threshold)

    return (by_platform_data)

class _markers():
  def __init__(self, db, filedata_coll, sample_coll, gwasdb, probidx=1):
    self.db = db
    self.gwasdb = gwasdb
    self.markers = db.markers
    self.calls = ["0/0", "0/1", "1/1", "Missing"]
    self.icalls = [0, 1, 2, -9]
    self.filedata_coll = filedata_coll
    self.sample_coll = sample_coll
    self.probidx = probidx
    self.vcfr = VCFrecord()

  def get_probidx(self):
    return self.probidx

  def get_marker_data(self, markerid, assaytype):
    """
    Get the data for a genetic marker / assay platform combination 
    """
    query = {}
    query['rsid'] = markerid
    query['assaytype'] = assaytype

    try:
      doc = self.markers.find_one(query)
    except:
      print "Unexpected error:", sys.exc_info()[0]
  
    # can return 'None' if query fails
    return doc

  def get_marker_data_multi(self, markerid):
    """
    Get the data for a genetic marker (DBSNP rs number or chrn:pos:I|D format)
    """
    query = {}
    query['rsid'] = markerid
    docs = []

    try:
      cursor = self.markers.find(query)
    except:
      print "Unexpected error:", sys.exc_info()[0]
  
    for doc in cursor:
      doc["samplecount"] = self.sample_coll.get_count(doc["assaytype"])
      docs.append(doc)
    # can return [] if query fails
    return docs

  def get_marker_data_by_range(self, chr, start, end):
    """
    Get the data for genetic markerwithin a range 
    """
    docs = []
    msg = ""
    start_pos = int(start)
    end_pos = int(end)

    # Some basic sanity checking
    if (end_pos - start_pos) > 250000:
      msg = "Range is too great should be 250Kb or less [%d]" % (end_pos - start_pos)
      return (docs, msg)
    if (end_pos - start_pos) < 0:
      msg = "Start pos is greater than End pos" 
      return (docs, msg)

    query = {}
    query['chromosome'] = chr = "%.2d" % (int(chr))
    query['position'] = {}
    query['position']['$gte'] = start_pos
    query['position']['$lte'] = end_pos

    print "RANGE QUERY", query

    try:
      cursor = self.markers.find(query)
    except:
      msg = "Unexpected error:" + sys.exc_info()[0]
  
    for doc in cursor:
      if len(doc["alleleA"]) > 10:
        doc["alleleA"] = doc["alleleA"][0:10] + " ..."
      if len(doc["alleleB"]) > 10:
        doc["alleleB"] = doc["alleleB"][0:10] + " ..."
      doc["samplecount"] = self.sample_coll.get_count(doc["assaytype"])
      docs.append(doc)
    # can return [] if query fails
    if len(docs) == 0:
      msg = "Nothing found in range"
    return (docs, msg)

  def get_marker_totals(self):
    """ Use agg framework to get totals by CHR
    """
    chr_totals = []
    curs = self.db.markers.aggregate([
      {"$group": {"_id":{"chr":"$chromosome", "at":"$assaytype"}, "mkrsPerChrom": { "$sum":1}}}, 
      {"$sort": {"_id":1}}
    ])
    for doc in curs['result']:
      #print doc['_id'], doc['mkrsPerChrom']
      chr_totals.append((doc['_id'], doc['mkrsPerChrom']))
    return chr_totals
  
  def get_raw_marker_values(self, filepath, variantid, chr, posn):
    """
    Access a vcf file to extract marker data
    """
    tabixFile = pysam.Tabixfile(filepath)
    if int(chr) > 22:
      chr = "NA"
    else:
      chr = str(chr)

    rec = []
    rtn_rec = ""

    try:
      records = tabixFile.fetch(chr, posn - 1, posn)
    except ValueError:
      chr = chr[1:]
      records = tabixFile.fetch(chr, posn - 1, posn)

    for record in records:
      data = self.vcfr.get_data_array(record)
      dvarid = self.vcfr.get_var_id_from_array(data)
      dposn = self.vcfr.get_posn_from_array(data)
      #print "%s-%s, %d-%d" % (variantid, dvarid, posn, int(dposn))
      if (dvarid == variantid) and (int(dposn) == posn):
        rec = data
        rtn_rec = record

    return rec, rtn_rec

  def get_genotype_probs(self, sample_values, threshold, has_GP=True):
    """
    Summarise marker_values based on probabilities
    TODO: deal with user-supplied threshold
    TODO: what to do when has_GP is false
    """
    geno_count = 0
    sample_count = 0
    genotype_counts = {}
    for sample_value in sample_values:
      sample_count += 1
      #print sample_value
      genoValues = sample_value.split(':')
      probVals = genoValues[self.get_probidx()].split(',')
      (key, ccode, maxprob) = self.get_call(probVals, threshold)
      genotype_counts[key] = genotype_counts.get(key, 0) + 1
    
    #print "SAMPLE COUNT", sample_count
    #gc_count_str = []
    hom1_ct = 0
    hom2_ct = 0
    het_ct = 0
    if "0/0" in genotype_counts:
      hom1_ct = genotype_counts["0/0"]
      geno_count += hom1_ct
    if "0/1" in genotype_counts:
      het_ct = genotype_counts["0/1"]
      geno_count += het_ct
    if "1/1" in genotype_counts:
      hom2_ct = genotype_counts["1/1"]
      geno_count += hom2_ct
    #print "allele counts:", hom1_ct, het_ct, hom2_ct, sample_count
    mafr = maf(het_ct, hom1_ct, hom2_ct, geno_count)
    #print "mafr:", mafr
    AlleleAfr = af(het_ct, hom1_ct, hom2_ct, geno_count)
    #print "afr:", AlleleAfr
    AlleleBfr = af(het_ct, hom2_ct, hom1_ct, geno_count)
    #print "bfr:", AlleleBfr
    p_hwe = HWE_exact(het_ct, hom1_ct, hom2_ct, geno_count)
    #for gt in genotype_counts:
    #  gc_count_str.append(gt + ": " + str(genotype_counts[gt])) 
    return (genotype_counts, sample_count, geno_count, mafr, AlleleAfr, AlleleBfr, p_hwe)

  def get_call(self, probs, threshold):
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
  
    return (self.calls[max_idx], self.icalls[max_idx], max_prob)

  def get_geno_data(self, rsid, sample_id, assaytype_list_posns):
    first_sample_idx = 9
    geno_values = {}
    docs = self.get_marker_data_multi(rsid)

    for doc in docs:
      # first get filepath
      fpath = self.filedata_coll.get_filepath(doc["assaytype"], doc["chromosome"])
      # this is not ideal - need to get on top of this chromosome id thing
      chr = "%.2d" % int(doc["chromosome"])
      rec, fullrec = self.get_raw_marker_values(fpath, doc["rsid"], chr, doc['position'])

      if doc["assaytype"] in assaytype_list_posns:
        prfx, genodata = self.vcfr.get_prfx_sfx_from_array(rec)
        print "list posn", rsid, sample_id, assaytype_list_posns[doc["assaytype"]]
        geno_values[sample_id + "_" + doc["assaytype"]] = genodata[assaytype_list_posns[doc["assaytype"]]]
    return (geno_values)


class _filedata():
  def __init__(self, db):
    self.db = db
    self.filedata = db.filedata

  def get_filepath(self, assaytype, chr):
    """
    Get the file to open for the genotype part of the query
    """
    query = {}
    query["assaytype"] = assaytype

    try:
      doc = self.filedata.find_one(query)
    except:
      print "Unexpected error:", sys.exc_info()[0]
  
  # can return 'None' if query fails
    if doc == None:
      return None
  
    filepath = str(doc["filepath"])
    rtn_file = ""

    for file in doc["files"]:
      if chr == str(file["CHROM"]):
        filepath = filepath + "/" + str(file["filename"])

    return filepath

class _samples():
  def __init__(self, db):
    self.db = db
    self.samples = db.samples
    self.sampleLists = {}
    self.sample_counts = {} # by assaytype

  def get_count(self, assaytype):
    samples = self.samples
    #print assaytype
    if assaytype not in self.sample_counts:
      try:
        query = {}
        query['assaytype'] = assaytype
        self.sample_counts[assaytype] = samples.find(query).count()
      except:
        print "Unexpected ERROR:", sys.exc_info()
    return self.sample_counts[assaytype]
  
  def get_samples(self, assaytype):
    """
    Get list of sampleids for assaytype in index order (that is, the order found in the .gen and .vcf files)
    """
    if assaytype in self.sampleLists:
      return self.sampleLists[assaytype]

    self.sampleLists[assaytype] = []
    try:
      query = {}
      query['assaytype'] = assaytype
      cursor = self.samples.find(query)
      cursor = cursor.sort([('list_posn',pymongo.ASCENDING)])
    except:
      print "Unexpected error:", sys.exc_info()[0]

    for doc in cursor:
      self.sampleLists[assaytype].append(doc["sample_id"])

    return self.sampleLists[assaytype]

  def get_all_samples(self):
    """
    Get list of all sampleids regardless of assaytype
    """
    samplelist = []
    try:
      query = {}
      cursor = self.samples.find(query)
    except:
      print "Unexpected error:", sys.exc_info()[0]

    for doc in cursor:
      samplelist.append(doc["sample_id"])

    return samplelist

  def get_sampledata(self, mprochid):
    """
    Get the list of samples collection docs for an mprochi
    """
    cursor=[]
    try:
      query = {}
      query["sample_id"] = mprochid
      cursor = self.samples.find(query)
    except:
      print "Unexpected error:", sys.exc_info()[0]

    docs = []
    for doc in cursor:
      docs.append(doc)
    return docs

class _prochi_map():
  def __init__(self, db, projpref, get_anochi):
    self.db = db
    self.prochi_map = db.prochi_map
    self.mprochiToProchi = {}
    self.plateIdToProchi = {}
    self.mprochiToAnochi = {}
    self.projpref = projpref
    self.get_anochi = get_anochi
    self.load_prochi_maps(projpref)

  def load_prochi_maps(self,projpref):
    """
    Load the two lookups Master ProCHI vs Prochi and plateId vs Prochi
    """
    query = {}
    cursor = self.prochi_map.find(query)

    for doc in cursor:
      if 'MProchi' in doc:
        if doc['Prochi'].startswith(projpref):
          self.mprochiToProchi[doc['MProchi']] = doc['Prochi']
      if 'Anochi' in doc:
        #print doc['MProchi'], doc['Anochi']
        self.mprochiToAnochi[doc['MProchi']] = doc['Anochi']
      if 'plateId' in doc:
        self.plateIdToProchi[doc['plateId']] = doc['Prochi']

  def get_anochi_or_prochi_from_mprochi(self, mprochi):
    """
    Get the either anochi or prochi
    """
    #print "Anochi logical = ", self.get_anochi, mprochi
    if self.get_anochi == True:
      return(self.get_anochi_from_mprochi(mprochi))
    else:
      return(self.get_prochi_from_mprochi(mprochi))

  def get_prochi_from_mprochi(self, mprochi):
    """
    Get the self.mprochiToProchi value for the supplied arg
    """
    try:
      return self.mprochiToProchi[mprochi]
    except:
      return None

  def get_anochi_from_mprochi(self, mprochi):
    """
    Get the self.mprochiToAnochi value for the supplied arg
    """
    try:
      print "Try", mprochi
      return self.mprochiToAnochi[mprochi]
    except:
      print "Except", mprochi
      return None

  def get_prochi_from_plateid(self, plateid):
    """
    Get the self.plateIdToProchi value for the supplied arg
    """
    try:
      return self.plateIdToProchi[plateid]
    except:
      return None

