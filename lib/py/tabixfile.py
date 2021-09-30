import logging
import pysam

class Tabixfile():
  def __init__(self):
    self.filepath = None
    self.tabixfiles = {}
    self.delim = "\t"

  def set_tabix_file(self, filepath):
    """
    Set the tabix file for multiple reads
    """
    #print "set_tabix_file", filepath
    if filepath not in self.tabixfiles:
      logging.info("Add tabix filepath - %s" % (filepath))
      self.tabixfiles[filepath] = pysam.Tabixfile(filepath)

  def get_file_data(self, filepath, chrom, start, end):
    """
    Access a vcf file to extract data, by chromosome, start pos, end pos
    """
    #print chrom, start, end
    res = []

    try:
      for record in self.tabixfiles[filepath].fetch(str(chrom), int(start) - 1, int(end)):
        res.append(record)
    except:
      logging.info("Tabix read error for - %s" % (self.filepath))

    return res

  def get_data_array(self, rec):
    return rec.split(self.delim)
