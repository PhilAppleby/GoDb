from tabixfile import Tabixfile

class Dbsnpfile(Tabixfile):

  def __init__(self):
    Tabixfile.__init__(self)
    self.rsid_idx = 2

  def get_dbsnp_file_record(self, filepath, chromosome, posn):
    recs = self.get_file_data(filepath, chromosome, posn, posn)
    for i, record in enumerate(recs):
      array = self.get_data_array(record)
      if int(array[1]) == int(posn):
        return record
    return None

  def get_rsid(self, rec):
    return self.get_data_array(rec)[self.rsid_idx]
