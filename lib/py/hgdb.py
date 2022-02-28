# Methods to access the human reference genome database (MongoDb)
import time
import pymongo

class Hgdb():
  def __init__(self, db, genome_str_len):
    self.db = db
    self.genome = db.genome
    self.genome_data = [] # store genome data for flushing to db
    self.genome_data_max_buffer = 1000
    self.genome_str_len = genome_str_len # shoould get this from the db
    self.add_count = 0
    self.first_genotype_idx = 5
    self.genotype_grp_size = 3

  def add_genome_data(self, chromosome, start_posn, genome_string, start_time):
    if len(self.genome_data) >= self.genome_data_max_buffer:
      self.flush_genome_buffer(start_time)
    
    doc = {}
    doc["chromosome"] = chromosome
    doc["start_posn"] = start_posn
    doc["genome_string"] = genome_string
    self.genome_data.append(doc)
    self.add_count += 1
  
  def flush_genome_buffer(self, start_time):
    self.genome.insert(self.genome_data)
    self.genome_data = [] # reset for next time
    #print ".", time.time() - start_time, "seconds", self.add_count

  def get_single_position(self, chromosome, offset):
    query = {}
    query["chromosome"] = chromosome
    query["start_posn"] = ((offset // self.genome_str_len) * self.genome_str_len) + 1
    #print query["start_posn"]
    index = offset - query["start_posn"]
    #print index

    try:
      doc = self.genome.find_one(query)
    except:
      #print "Unexpected error:", sys.exc_info()[0]
      return None
       
    # can return 'None' if query fails
    if doc == None:
      return None
     
    return doc["genome_string"][int(index)]

  def get_genotypes(self, full_record):
    """
    The argument is a split gen record
    """
    return (full_record[self.first_genotype_idx:])

  def parse_gen_prefix(self, full_record):
    """
    The argument is a split gen record
    """
    return (full_record[0], full_record[1], full_record[2], full_record[3], full_record[4])

  def flip(self, genotypes):
    """
    The argument is a genotype array 
    """
    for i in xrange(0, len(genotypes), self.genotype_grp_size):
      geno = genotypes[i]
      genotypes[i] = genotypes[i+2]
      genotypes[i+2] = geno
    return (genotypes)


