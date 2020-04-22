import pymongo
import pysam
from dbconfig import DBSERVER
from dbconfig import DBNAME
from dbconfig import VARIANTS
import re
import os, sys
from vcfrecord import VCFrecord
#
# Collection of methods to assist with GoDb management
# Loading methods for variants, samples, 
# filedata, filepaths and genemap 
# (filedata. filepaths are alternative versions for file location data)
#
class GoDb():
  def __init__(self):
    try:
      connection = pymongo.MongoClient(DBSERVER)
      db = eval("connection." + DBNAME)
      variants = eval("db." + VARIANTS)
    except:
      raise Exception("Unexpected error connecting to %s @ %s" % (DBNAME, DBSERVER))

    self.db = db
    self.variants = variants
    self.samples = db.samples
    self.filedata = db.filedata
    self.filepaths = db.filepaths
    self.genemap = db.genemap
    self.dbname = DBNAME
    self.variantbuff = []
    self.samplebuff = []
    self.filebuff = []
    self.genemapbuff = []
    self.int_fields = ["position"]
    self.flt_fields = ["all_maf", "info", "cohort_1_hwe"]
    self.p = re.compile('\d+')

  def get_dbname(self):
    return self.dbname
  
  # methods relating to the genemap collection
  
  def process_genemap_detail(self, record):
    """Process  ucsc refFlat file records
       Set up a json document and add it to the 
      genemap buffer
      format: 0;genename, 1;name, 2;chrom (in chr<N> format), 3;strand,
      4;txStart, 5;txEnd, 6:cdsStart, 7;cdsEnd, 8:exonCount, 9:exonStarts, 10:exonEnds
    """
    data = record.split()
    doc = {}
    exonStarts = []
    exonEnds = []

    doc["genename"] = data[0]
    doc["name"] = data[1]
    doc["chrom"] = data[2][3:]
    doc["strand"] = data[3]
    doc["txStart"] = int(data[4])
    doc["txEnd"] = int(data[5])
    doc["cdsStart"] = int(data[6])
    doc["cdsEnd"] = int(data[7])
    doc["exonCount"] = int(data[8])
    exStarts = data[9].split(",")
    for elem in exStarts:
      try:
        exonStarts.append(int(elem))
      except:
        pass
    doc["exonStarts"] = exonStarts
    exEnds = data[10].split(",")
    for elem in exEnds:
      try:
        exonEnds.append(int(elem))
      except:
        pass
    doc["exonEnds"] = exonEnds

    self.genemapbuff.append(doc)

  def flush_genemap_buff(self):
    try:
      if (len(self.genemapbuff) > 0):
        self.genemap.insert(self.genemapbuff)
    except:
      print "Unexpected error writing to genemap collection:", sys.exc_info()[0]
      sys.exit()

    self.genemapbuff = []

  def get_genemap_len(self):
    return len(self.genemapbuff)

  def get_one_genemap(self, genename):
    query = {}
    query["genename"] = genename

    try:
      doc = self.genemap.find_one(query)
    except:
      print "Unexpected error (get_one_genemap):", sys.exc_info()[0]
      sys.exit()
       
    # can return 'None' if query fails
    return doc

  def get_genemap_data_by_buffered_range(self, chrom, posn, upstream, downstream):
    query = {}
    docs = []
    # This is effectively widening the range checked by shifting the
    # genomic postion being testes
    query['txStart'] = {}
    query['txEnd'] = {}
    query['chrom'] = chrom
    query['txStart']['$lte'] = posn + upstream
    query['txEnd']['$gte'] = posn - downstream

    try:
      curs = self.genemap.find(query)
      curs = curs.sort([('genename',pymongo.ASCENDING), ('txStart',pymongo.ASCENDING), ('txEnd',pymongo.ASCENDING)])
    except:
      print "Unexpected error (get_genemap_data_by_buffered_range):", sys.exc_info()[0]
      sys.exit()
    for doc in curs:
      #print "%s,%s,%d,%d" % (doc["genename"], doc["chrom"], doc["txStart"], doc["txEnd"])
      docs.append(doc)
    # can return 'None' if query fails
    return docs

  # methods relating to the variants collection

  def process_variant_detail(self, hdr, record, assaytype):
    """Process info file variant detail records
       Set up a json document and add it to the 
      variant buffer
    """
    doc = {}
    doc["assaytype"] = assaytype
    for elem_num, field in enumerate(hdr):
      if field in self.int_fields:
        doc[field] = int(record[elem_num])
      elif field in self.flt_fields:
        doc[field] = float(record[elem_num])
      elif record[elem_num] == "---":
        doc["imputed"] = 1
      else:
        doc[field] = record[elem_num]

    self.variantbuff.append(doc)

  def process_variant_detail_vcf(self, record, assaytype):
    """Process info file variant detail records
       Set up a json-stype document and add it to the 
      variant buffer
    """
    doc = {}
    doc["assaytype"] = assaytype
    vcfr = VCFrecord(record)
    prfx, sfx = vcfr.get_prfx_sfx() 
    doc["rsid"] = vcfr.get_varid()  
    # always store chromosome as a 2-digit string
    doc["chromosome"] = "%.2d" % (int(vcfr.get_chr()))
    alleleA, alleleB = vcfr.get_alleles()
    doc["alleleA"] = alleleA
    doc["alleleB"] = alleleB
    doc["position"] = vcfr.get_posn_as_int()
    try:
      doc["ref_maf"] = float(vcfr.get_info_value("RefPanelAF"))  
    except:
      pass
    try:
      doc["info"] = float(vcfr.get_info_value("INFO"))  
    except:
      doc["info"] = 1.0

    self.variantbuff.append(doc)

  def flush_variant_buff(self):
    try:
      if (len(self.variantbuff) > 0):
        self.variants.insert(self.variantbuff)
    except:
      print "Unexpected error writing to variant collection:", sys.exc_info()[0]
      sys.exit()

    self.variantbuff = []

  def get_variants_len(self):
    return len(self.variantbuff)

  def get_one_variant_for_assaytype(self, rsid, assaytype):
    """
    Find one variant "variant" record, for an rsid, assaytype
    """
    query = {}
    query["rsid"] = rsid
    query["assaytype"] = assaytype

    try:
      doc = self.variants.find_one(query)
    except:
      print "Unexpected error:", sys.exc_info()[0]
      sys.exit()
       
    # can return 'None' if query fails
    return doc

  def get_one_variant(self, rsid):
    query = {}
    query["rsid"] = rsid

    try:
      doc = self.variants.find_one(query)
    except:
      print "Unexpected error:", sys.exc_info()[0]
      sys.exit()
       
    # can return 'None' if query fails
    return doc

  def get_multiple_variants(self, rsid):
    query = {}
    query["rsid"] = rsid

    try:
      curs = self.variants.find(query)
    except:
      print "Error:", sys.exc_info()[0]
      sys.exit()
       
    # can return 'None' if query fails
    return curs

  def get_variant_data_by_range(self, chromosome, start, end):
    """ 
    Get the data for variants within a genomic range 
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
    query['chromosome'] = chromosome = "%.2d" % (int(chromosome))
    query['position'] = {}
    query['position']['$gte'] = start_pos
    query['position']['$lte'] = end_pos

    try:
      cursor = self.variants.find(query)
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

  # methods relating to the samples collection

  def process_sample_detail(self, sample_id, idx, assaytype):
    """Process sample date (called once per sample
       Set up a json-type document and add it to the 
      sample buffer
    """
    doc = {}
    doc["sample_id"] = sample_id
    doc["assaytype"] = assaytype
    doc["list_posn"] = idx

    self.samplebuff.append(doc)

  def flush_sample_buff(self):
    try:
      if (len(self.samplebuff) > 0):
        self.samples.insert(self.samplebuff)
    except:
      print "Unexpected error writing to samples collection:", sys.exc_info()[0]
      sys.exit()

    self.samplebuff = []

  def get_samples_len(self):
    return len(self.samplebuff)
  
  def get_samples(self, assaytype):
    """
    Get the list of samples for an assaytype in list_posn order
    """
    sample_list = []
    try:
      query = {}
      query['assaytype'] = assaytype
      cursor = self.samples.find(query)
      cursor = cursor.sort([('list_posn',pymongo.ASCENDING)])
    except:
      print "Unexpected error:", sys.exc_info()[0]
      sys.exit()

    for doc in cursor:
      sample_list.append(doc["sample_id"])

    return sample_list
  
  def get_sample_count(self, assaytype):
    """
    Get the count of samples for an assaytype 
    """
    count = 0
    try:
      query = {}
      query['assaytype'] = assaytype
      count = self.samples.find(query).count()
    except:
      print "Unexpected ERROR:", sys.exc_info()
      sys.exit()

    return count

  # methods relating to the filepaths collection

  def add_filepath_detail(self, assaytype, fprefix, fsuffix, filelist):
    """Process filepath date (called for each file set).
    """
    doc = {}
    doc["assaytype"] = assaytype
    doc["fpath_prefix"] = fprefix
    doc["fpath_suffix"] = fsuffix
    doc["filepath"] = fprefix + "/" + fsuffix

    files = []
    for filename in filelist:
      m = self.p.search(filename)
      if m != None:
        filedoc = {}
        filedoc["CHROM"] = m.group()
        filedoc["filename"] = filename
        files.append(filedoc)

    doc["files"] = files

    try:
      self.filepaths.insert(doc)
    except:
      print "Unexpected error writing to filepaths collection:", sys.exc_info()[0]
      sys.exit()

  def get_filepath(self, assaytype, chromosome):
    query = {}
    query["assaytype"] = assaytype

    try:
      doc = self.filepaths.find_one(query)
    except:
      print "Unexpected error:", sys.exc_info()[0]

    if doc == None:
      return None    

    filepath = str(doc["filepath"])

    for filedoc in doc["files"]:
      if chromosome == str(filedoc["CHROM"]):
        filepath = filepath + "/" + str(filedoc["filename"])

    return filepath

  def get_full_filepath(self, assaytype, chromosome, prfx):
    """ 
    """
    query = {}
    query["assaytype"] = assaytype
    
    try:
      doc = self.filepaths.find_one(query)
    except:
      print "Unexpected error:", sys.exc_info()[0]
    
    # can return 'None' if query fails
    if doc == None:
      return None
    
    filepath_suff = str(doc["fpath_suffix"])
    filepath = prfx + "/" + filepath_suff
    rtn_file = ""

    for file in doc["files"]:
      if chromosome == str(file["CHROM"]):
        filepath = filepath + "/" + str(file["filename"])

    return filepath



  # methods relating to tabix files
  def get_variant_file_data(self, filepath, chromosome, posn):
    """ 
    Access a vcf file to extract variant data
    """
    tabixFile = pysam.Tabixfile(filepath)
    if int(chromosome) > 22: 
      chromosome = "NA"
    else:
      chromosome = str(chromosome)

    rtn_rec = None

    try:
      records = tabixFile.fetch(chromosome, posn - 1, posn)
    except ValueError:
      chromosome = chromosome[1:]
      records = tabixFile.fetch(chromosome, posn - 1, posn)

    for record in records:
      rtn_rec = record

    return rtn_rec

