import logging
import sys
import re
from hwehelper import Hwehelper
from mafhelper import Mafhelper
from godb import GoDb
from vcfrecord import VCFrecord
from zipfile import ZipFile
from zipfile import ZIP_DEFLATED
from config import RSLISTMAX

class DataStore():
  def __init__(self):
    self.godb = GoDb()
    self.filepaths_coll = _filepaths(self.godb)
    self.sam_coll = _samples(self.godb)
    self.var_coll = _variants(self.filepaths_coll, self.sam_coll, self.godb)
    self.variant_totals = []
    self.sample_count = -1
    self.call_rates = {}
    self.maxrslist = 200

  def get_db_name(self):
    return self.godb.get_dbname()

  def make_selection_key(self, varid, assaytype):
    return varid + "_" + assaytype

  def get_variant_data_for_file(self, filepath, threshold):
    msg = ""
    pattern = re.compile('[\W_]+')
     
    try:
      f = open(filepath, "r")
    except IOError as e:
      msg = filepath + ":" + e.strerror
      return([],msg)
    count = 0
    
    rslist = []
    for line in f:
      count += 1
      if count > self.maxrslist:
        msg = "line count for %s gt the limit (%d)" % (filepath, self.maxrslist)  
        return ([], msg)
      line = line.strip()
      pattern.sub('', line)
      if line.startswith("rs"):
        elems = line.split()
        rslist.append(elems[0])
      else:
        msg += "Removed bad line at %d, " % (count)

    f.close()
    variant_data = []
    for rsid in rslist:
      (variant_docs, tmsg) = self.get_variant_summary_probs(rsid, threshold)
      if (tmsg != ""):
        msg+=tmsg
      for doc in variant_docs:
        variant_data.append(doc)
    return variant_data, msg

  def get_variant_data_by_range(self, chromosome, start, end, threshold=0.9):
    rslist = []
    variant_data = []
    msg = ""
    return(self.var_coll.get_variant_data_by_range(chromosome, start, end))

  def get_range_data(self, chromosome, start, end, threshold, download_list):
    (docs, msg) = self.var_coll.get_variant_data_by_range(chromosome, start, end)
    if len(docs) == 0:
      return([], [], msg)

    rsdict={}
    rslist=[]
    for doc in docs:
      rsdict[doc["rsid"]] = 1

    for rsid in rsdict:
      rslist.append(rsid)

    return(self.get_rslist_data(rslist, threshold, download_list))

  def build_csv_data(self, rslist, sampleDict, assaytypes, threshold, atpidx):
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
          geno_data = self.resolve_geno(sampleDict[samp][rsid], rsid, samp, threshold, atpidx)
          dataVals = geno_data[2].split(':')
          #print "build_csv_data:", geno_data
          probVals = dataVals[atpidx[geno_data[0]]].split(',')
          (probcall, intcall, outprob) = self.var_coll.get_call(probVals, threshold)
          output_lines['combined'][idxoffset] = str(intcall)
          output_lines['combined'][idxoffset+1] = str(outprob)
          output_lines['combined'][idxoffset+2] = geno_data[0]
          output_lines['combined'][idxoffset+3] = geno_data[4]
          for geno_data in sampleDict[samp][rsid]:
            dataVals = geno_data[2].split(':')
            probVals = dataVals[atpidx[geno_data[0]]].split(',')
            (probcall, intcall, outprob) = self.var_coll.get_call(probVals, threshold)
            output_lines[geno_data[0]][idxoffset] = str(intcall)
            output_lines[geno_data[0]][idxoffset+1] = str(outprob)
            output_lines[geno_data[0]][idxoffset+2] = geno_data[0]
            output_lines[geno_data[0]][idxoffset+3] = geno_data[4]
            filled_output_lines[geno_data[0]] = True
      
      for assaytype in assaytypes:
        if assaytype in filled_output_lines:
          by_platform_data[assaytype].append(samp + "," + ",".join(output_lines[assaytype]))
    return by_platform_data
  
  def resolve_geno(self, genlist, rsid, samp, threshold, atpidx):
    maxprob = 0.0
    maxidx = -1
    if len(genlist) == 1:
      return genlist[0]
    elif len(genlist) > 1:
      for idx, gendata in enumerate(genlist):
        try:
          dataVals = gendata[2].split(':')
          probVals = dataVals[atpidx[gendata[0]]].split(',')
        except IndexError:
          sys.exit()
        (probcall, intcall, outprob) = self.var_coll.get_call(probVals, threshold)
        if outprob > maxprob:
          maxprob = outprob
          maxidx = idx
    if maxprob > 0.0:
      return genlist[maxidx]
    return genlist[0]

  def get_variant_summary_probs(self, rsid, threshold):
    variant_array = []
    msg = ""
    docs = self.var_coll.get_variant_data_multi(rsid)
    for doc in docs:
      # always force chromosome to 2 digits
      chromosome = "%.2d" % int(doc["chromosome"])
      fpath = self.filepaths_coll.get_filepath(doc["assaytype"], chromosome)
      fullrec = self.var_coll.get_raw_variant_values(fpath, chromosome, doc['position'])
      vcfr = VCFrecord(fullrec)
      prfx, sfx = vcfr.get_prfx_sfx()
      probidx = vcfr.get_probidx()
      (gc_count_dict, sample_count, geno_count, maf, alleleAf, alleleBf, p_hwe) = self.var_coll.get_genotype_probs(sfx, threshold, probidx)
      doc['selected'] = 1
      doc['a_af'] = alleleAf
      doc['b_af'] = alleleBf
      doc['hwe_p'] = p_hwe
      doc['Missing'] = 0
      if 'Missing' in gc_count_dict:
        doc['Missing'] = gc_count_dict['Missing']
      variant_array.append(doc)
    if len(variant_array) == 0:
      msg = "Variant NOT FOUND - %s, " % (rsid)
    return (variant_array, msg)

  def get_variant_data(self, variantid, assaytype):
    return self.var_coll.get_variant_data(variantid, assaytype)

  def get_sample(self, sampleid):
    return self.sam_coll.get_sample(sampleid)

  def make_zipfile(self, sample_return_data, snp_return_data, uploadDir, zipfilename):
    """
    moved here from views.py
    """
    ares = {}
    for assaytype in sample_return_data:
      ares[assaytype] = '\n'.join(sample_return_data[assaytype]) + '\n'
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
    assaytypelist = []
    probidxlist = []
    rslist = []
    assaytypes = {}
    Afreq = {}
    Bfreq = {}

    data_count = 0
    impDict = {}
    for rsid in input_rslist:
      docs = self.var_coll.get_variant_data_multi(rsid)
      if len(docs) > 0:
        rslist.append(rsid)
      # handling SNPs on multiple platforms
      for doc in docs:
        # always force chromosome to 2 digits
        chromosome = "%.2d" % int(doc["chromosome"])
        # first get filepath
        fpath = self.filepaths_coll.get_filepath(doc["assaytype"], chromosome)
        # get raw variant data
        fullrec = self.var_coll.get_raw_variant_values(fpath, chromosome, doc['position'])
        geno_count = 0
        sample_count = 0
        hwep = 0.0
        vcfr = VCFrecord(fullrec)
        prfx, sfx = vcfr.get_prfx_sfx()
        probidx = vcfr.get_probidx()
        (gc_count_dict, sample_count, geno_count, maf, alleleAf, alleleBf, p_hwe) = self.var_coll.get_genotype_probs(sfx, threshold, probidx)
        Afreq[doc["rsid"] + "_" + doc["assaytype"]] = alleleAf
        Bfreq[doc["rsid"] + "_" + doc["assaytype"]] = alleleBf
        data_count += 1
        hwep = float(p_hwe) 
          
        assaytypelist.append(doc["assaytype"])
        data.append(vcfr)
        if doc["assaytype"] not in assaytypes:
          assaytypes[doc["assaytype"]] = 1

        imputed = 0
        if "imputed" in doc:
          imputed = 1
        if "info" in doc:
          if doc["info"] != 1.0:
            imputed = 1
        impDict[doc["rsid"] + "_" + doc["assaytype"]] = imputed
        snpdata += "%s,%s,%s,%d,%s,%s,%s,%s,%s,%d,%.5f,%.5f,%.5f\n" % (doc["rsid"], doc["assaytype"], doc["chromosome"], doc["position"], doc["alleleA"], alleleAf, doc["alleleB"], alleleBf, maf, imputed, float(geno_count) / sample_count, hwep, doc["info"])

    pdata = self.get_sample_values(assaytypelist, data, data_count, rslist, impDict, assaytypes, Afreq, Bfreq, threshold)
    return(pdata, snpdata, msg)

  def get_sample_values(self, assaytypelist, records, numrecs, rslist, impDict, 
      assaytypes, Afreq, Bfreq, threshold):
    """Process vcf records
      NOTE: impDict keyed on a composite of rsid_platform 
       """
    samplesByAt = {}
    for assaytype in assaytypes:
      samplesByAt[assaytype] = self.sam_coll.get_samples(assaytype)
    # A dict of dicts of tables
    sampleDict = {}
    dupDict = {}
    dupcount = 0
    for assaytype in samplesByAt:
      for samp in samplesByAt[assaytype]:
        if samp not in sampleDict:
          sampleDict[samp] = {}
          for rsid in rslist:
            sampleDict[samp][rsid] = []

    count = 0
    rscount = 0
    atpidx = {}
    values_totals = 0
    for idx, vcfr in enumerate(records):
      assaytype = assaytypelist[idx]
      chromosome = vcfr.get_chr()
      pos = vcfr.get_posn()
      rsid = vcfr.get_varid()
      alleleA, alleleB = vcfr.get_alleles()
      prfx, sfx = vcfr.get_prfx_sfx()
      atpidx[assaytype] = vcfr.get_probidx()

      if assaytype not in samplesByAt:
        samplesByAt[assaytype] = self.sam_coll.get_samples(assaytype)

      for idx, elem in enumerate(sfx):
        if samplesByAt[assaytype][idx] in sampleDict:
          sampleId = samplesByAt[assaytype][idx]
          values_totals += 1
          sampleDict[sampleId][rsid].append([assaytype,impDict[rsid  
            + "_" + assaytype],sfx[idx],alleleA,alleleB, 
            Afreq[rsid  + "_" + assaytype], Bfreq[rsid  + "_" + assaytype]])
      rscount += 1 

    by_assaytype_data = self.build_csv_data(rslist, sampleDict, assaytypes, threshold, atpidx)

    return (by_assaytype_data)

class _variants():
  def __init__(self, filepaths_coll, sample_coll, godb):
    self.godb = godb
    self.calls = ["0/0", "0/1", "1/1", "Missing"]
    self.icalls = [0, 1, 2, -9]
    self.filepaths_coll = filepaths_coll
    self.sample_coll = sample_coll
    self.hweh = Hwehelper()
    self.mafh = Mafhelper()

  def get_variant_data(self, variantid, assaytype):
    """
    Get the data for a genetic variant / assay platform combination 
    """
    doc = self.godb.get_one_variant_for_assaytype(assaytype)
  
    # can return 'None' if query fails
    return doc

  def get_variant_data_multi(self, variantid):
    """
    Get the data for a genetic variant (DBSNP rs number or chrn:pos:I|D format)
    """
    docs = []

    cursor = self.godb.get_multiple_variants(variantid)
  
    for doc in cursor:
      doc["samplecount"] = self.sample_coll.get_count(doc["assaytype"])
      docs.append(doc)
    # can return [] if query fails
    return docs

  def get_variant_data_by_range(self, chromosome, start, end):
    """
    Get the data for variants within a genomic range 
    """
    # delegate to godb
    return self.godb.get_variant_data_by_range(chromosome, start, end)

  def get_raw_variant_values(self, filepath, chromosome, posn):
    """
    Access a vcf file to extract variant data
    """
    return self.godb.get_variant_file_data(filepath, chromosome, posn)

  def get_genotype_probs(self, sample_values, threshold, probidx, has_GP=True):
    """
    Summarise variant_values based on probabilities
    TODO: deal with user-supplied threshold
    TODO: what to do when has_GP is false
    """
    geno_count = 0
    sample_count = 0
    genotype_counts = {}
    for sample_value in sample_values:
      sample_count += 1
      genoValues = sample_value.split(':')
      probVals = genoValues[probidx].split(',')
      (key, ccode, maxprob) = self.get_call(probVals, threshold)
      genotype_counts[key] = genotype_counts.get(key, 0) + 1
    
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
    mafr, ma = self.mafh.maf(het_ct, hom1_ct, "", hom2_ct, "", geno_count)
    AlleleAfr = self.mafh.af(het_ct, hom1_ct, hom2_ct, geno_count)
    AlleleBfr = self.mafh.af(het_ct, hom2_ct, hom1_ct, geno_count)
    p_hwe = self.hweh.HWE_exact(het_ct, hom1_ct, hom2_ct, geno_count)
    return (genotype_counts, sample_count, geno_count, mafr, AlleleAfr, AlleleBfr, p_hwe)

  def get_call(self, probs, threshold):
    max_prob = 0.0
    max_idx = 3

    for idx, prob in enumerate(probs):
      if float(prob) > max_prob:
        max_prob = float(prob)
        max_idx = idx

    if (threshold !=0.0):
      if max_prob < threshold:
        max_idx = 3
  
    return (self.calls[max_idx], self.icalls[max_idx], max_prob)

  def get_geno_data(self, rsid, sample_id, assaytype_list_posns):
    geno_values = {}
    docs = self.get_variant_data_multi(rsid)

    for doc in docs:
      # always force chromosome to 2 digits
      chromosome = "%.2d" % int(doc["chromosome"])
      fpath = self.filepaths_coll.get_filepath(doc["assaytype"], chromosome)
      fullrec = self.get_raw_variant_values(fpath, chromosome, doc['position'])

      if doc["assaytype"] in assaytype_list_posns:
        vcfr = VCFrecord(fullrec)
        prfx, genodata = vcfr.get_prfx_sfx()
        geno_values[sample_id + "_" + doc["assaytype"]] = genodata[assaytype_list_posns[doc["assaytype"]]]
    return (geno_values)


class _filepaths():
  def __init__(self, godb):
    self.godb = godb

  def get_filepath(self, assaytype, chromosome):
    """
    Get the file to open for the genotype part of the query
    """
    filepath = self.godb.get_filepath(assaytype, chromosome)
  # can return 'None' if query fails
    return filepath

class _samples():
  def __init__(self, godb):
    self.godb = godb
    self.sampleLists = {}
    self.sample_counts = {} # by assaytype

  def get_count(self, assaytype):
    if assaytype not in self.sample_counts:
      self.sample_counts[assaytype] = self.godb.get_sample_count(assaytype)

    return self.sample_counts[assaytype]
  
  def get_samples(self, assaytype):
    """
    Get list of sampleids for assaytype in index order (that is, the order found in the .gen and .vcf files)
    """
    if assaytype in self.sampleLists:
      return self.sampleLists[assaytype]

    self.sampleLists[assaytype] = self.godb.get_samples(assaytype)

    return self.sampleLists[assaytype]

