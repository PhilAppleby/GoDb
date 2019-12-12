import logging
from hwehelper import Hwehelper
from mafhelper import Mafhelper
from vcfrecord import VCFrecord

from multimerge import Multimerge
# Helper methods for merging multiple VCF files simultaneously
class Multifilemerge(Multimerge):
  def __init__(self, fh_list):
    Multimerge.__init__(self, len(fh_list), vcfr)
    self.fh_list = fh_list
    self.mafh = Mafhelper()
    self.prt_comments = []
    self._load_header_data()
    self._load_col_data()
    self.high_key = 999999999999999
    self.low_key = 0
    self.empty_key = -1
    self.rec_counts = [0] * self.numobjects
    self.write_count = 0
    self.allele_discord_count = 0
    self.gt1_match_count = 0
    self.not_concordant_tot = 0
    self.hweh = Hwehelper()
    self.col_hdr_pref = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL",  "FILTER",  "INFO",  "FORMAT"]
    self.at_header = "##FORMAT=<ID=AT,Number=1,Type=String,Description=\"Assay Type\">"

  def _load_header_data(self):
    """
    Not meant to be public - called at init time
    Load all the data contained in the header records.
    Includes the "##" information content and the "#" column headers
    """
    for i, fh in enumerate(self.fh_list):
      get_cmts = True
      while get_cmts == True:
        line = fh.readline().strip()
        if line.startswith("##"): 
          if i == 0:
            self.prt_comments.append(line)
        elif line.startswith("#"):
          vcfr = VCFRecord(line)
          prfx, sfx = vcfr.get_prfx_sfx()
          self.hdr_cols[i] = sfx
          get_cmts = False
        else:
          get_cmts = False

    return 

  def get_header_rows(self):
    return self.prt_comments

  def write_headers(self, fh_out):
    for cmt in self.prt_comments:
      fh_out.write(cmt + "\n")

    fh_out.write(self.at_header + "\n")

  def write_col_headers(self, fh_out):
    fh_out.write("\t".join(self.col_hdr_pref + self.get_combined_columns()) + "\n")

  def get_init_key_list(self):
    return [self.low_key] * self.numobjects

  def get_empty_buffers(self):
    return [[]] * self.numobjects
  
  def get_low_key_list(self, key_list):
    """ 
    The argument is a list of keys, which are ints
    """
    low_key = self.high_key

    for key in key_list:
      if key < low_key:
        low_key = key 

    rtn_keys = [self.empty_key] * self.numobjects

    low_key_count = 0
    for i, key in enumerate(key_list):
      if key == low_key:
        low_key_count += 1
        rtn_keys[i] = key 

    return rtn_keys, low_key_count

  def write_from_buffers(self, fh_out, key_list, prfx_list, recbuff_list, assaytype_list, atype_posns, chipval):
    """
    main rule is we only write from buffers corresponding to the minimum keys (we get the indexes for these by calling self.get_low_keys()
    """
    low_key_list, low_key_count = self.get_low_key_list(key_list)

    outbuff_list = [[]] * len(assaytype_list)

    prfx = ""

    for i, key in enumerate(low_key_list):
      if key != self.empty_key:
        outbuff_list[i] = recbuff_list[i] 
        prfx = prfx_list[i]

    concordant = True
    if low_key_count > 1:
      self.gt1_match_count += 1
      concordant = self.check_concordancies(prfx_list, outbuff_list, assaytype_list, chipval)

    if concordant == True:
      probidx = vcfr.get_probidx_from_array(prfx)
      comborec = self.get_combined_array(outbuff_list, assaytype_list, atype_posns, probidx)
      prfx[8] += ":AT"
      outrec = prfx + comborec
      self.write_count += 1
      return fh_out.write("\t".join(outrec) + "\n")
    else:
      self.not_concordant_tot += 1
      return 0

  def get_next_records(self, key_list, prfx_list, recbuff_list):
    """
    main rule is we read from the fh's corresponding to the min key list and replace the key_list, prfx and rec_buff elements accordingly.
    """
    low_key_list, low_key_count = self.get_low_key_list(key_list)

    for i, fh in enumerate(self.fh_list):
      if low_key_list[i] != self.empty_key:
        line = fh.readline().strip()
        if line != "": # testing for EOF
          self.rec_counts[i] += 1
          vcfr = VCFrecord(line)
          prfx, sfx = vcfr.get_prfx_sfx()
          maf, ma, cr = self.mafh.get_maf_and_cr(data, vcfr)           
          prfx_list[i] = prfx
          recbuff_list[i] = sfx
          key_list[i] = int(prfx[1])
        else:
          prfx_list[i] = []
          recbuff_list[i] = []
          key_list[i] = self.high_key
    #logging.info("rec_counts: %s, key_list: %s" % (str(self.rec_counts), str(key_list)))
    return key_list, prfx_list, recbuff_list

  def all_files_finished(self, key_list):

    for key in key_list:
      if key != self.high_key:
        return False

    return True

  def get_rec_counts(self):
    return self.rec_counts, self.write_count, self.chisq_count, self.allele_discord_count, self.gt1_match_count, self.not_concordant_tot


