import logging
from scipy.stats import chisquare # Used in hwe concordancy check
from vcfrecord import VCFrecord 
# Helper methods for merging data copies simultaneously
# Assumptions:
#   
#   Many of the instance variables are collections of 
#   size len(fh_list) apart from the combined ones which are 
#   "flattened"
# Loading of collections the responsibility of subclasses
#
class Multimerge():
  def __init__(self, numobjects):
    self.numobjects = numobjects
    self.file_positions = [dict() for x in range(self.numobjects)]
    self.combined_positions = {} 
    self.combined_columns = [] 
    self.hdr_cols = [[]] * self.numobjects
    self.assay_abbrev = {}
    self.assay_abbrev["affy"] = "A"
    self.assay_abbrev["illumina"] = "I"
    self.assay_abbrev["exome"] = "E"
    self.assay_abbrev["metabo"] = "M"
    self.assay_abbrev["affy1KG"] = "A1"
    self.assay_abbrev["illumina1KG"] = "I1"
    self.assay_abbrev["broad"] = "B"
    self.assay_abbrev["bigtest"] = "T"
    self.assay_abbrev["biggertest"] = "G"
    self.assay_expand = {}
    self.assay_expand["A"] = "affy"
    self.assay_expand["I"] = "illumina"
    self.assay_expand["E"] = "exome"
    self.assay_expand["M"] = "metabo"
    self.assay_expand["A1"] = "affy1KG"
    self.assay_expand["I1"] = "illumina1KG"
    self.assay_expand["B"] = "broad"
    self.assay_expand["T"] = "bigtest"
    self.assay_expand["G"] = "biggertest"
    self.geno_strings = ["0/0", "0/1", "1/1"]
    self.chisq_count = 0
    self.allele_discord_count = 0
    self.geno_overlap_count = 0

  def get_counts(self):
    return self.chisq_count, self.allele_discord_count, self.geno_overlap_count

  def _load_col_data(self):
    """ 
    Not meant to be public - called at init time
    Load the column list
    """
    z = 0 
    combo_temp = {}
    for i, colname_list in enumerate(self.hdr_cols):
      for y, colname in enumerate(colname_list):
        self.file_positions[i][y] = colname 
        if colname not in self.combined_positions:
          self.combined_positions[colname] = z 
          combo_temp[z] = colname 
          z += 1
    for key in sorted(combo_temp):
      self.combined_columns.append(combo_temp[key])
    return 

  def get_header_cols(self):
    return self.hdr_cols

  def get_combined_positions(self):
    return self.combined_positions

  def get_combined_columns(self):
    return self.combined_columns

  def get_combined_array(self, buffer_list, cr_list, assay_list, threshold=0.9):
    """
    For each list of data, for each element of list of data:
    1) Find the col header from the corresonding file_position element
    2) Use the col_header to find the combined postion
    3) Place the data_element in the combined postion *
    TODO - conflict resolution, what to do if a slot is already occupied
    TODO - CR check
    """
    #print "COMBO", self.combined_positions
    #
#print "ASSAY_LIST: %s" % (str(assay_list))
    assay_posns = {}

    for i, assaytype in enumerate(assay_list):
      assay_posns[i] = assaytype

    #print "ASSAY_POSNS: %s" % (str(assay_posns))

    combo_array = ["."] * len(self.combined_positions)
    #print "BUFFL", len(buffer_list)
    for i, vcf_record in enumerate(buffer_list): 
      if len(vcf_record) > 0:
        #print "asstp: %d, %s" % (i, assay_list[i])
        vcfr = VCFrecord(vcf_record)
        prfx, data_list = vcfr.get_prfx_sfx()
        probidx = vcfr.get_probidx()
        rsid = vcfr.get_varid()
        for j, dataelem in enumerate(data_list):
          if data_list[j] != ".":
            cpos = self.combined_positions[self.file_positions[i][j]]
            geno = self.call_geno_for_threshold(data_list[j], probidx, threshold) + ":" + self.assay_abbrev[assay_list[i]]
            if combo_array[cpos] != ".":
              self.geno_overlap_count += 1
              #print "OVERLAP %s:%s - %s vs %s" % (rsid, self.file_positions[i][j], combo_array[cpos], geno)
              geno = self.call_genotype(combo_array[cpos], geno, probidx)
            combo_array[cpos] = geno

    return combo_array


  def call_genotype(self, geno1, geno2, probidx=1):
    gen_vals1 = geno1.split(":")
    gen_vals2 = geno2.split(":")
    #if gen_vals1[0] == gen_vals2[0]:
    #  logging.info("Call geno EQUAL:%s vs %s at %d" % (geno1, geno2, cpos))
    #else:
    #  logging.info("Call geno DIFF:%s vs %s at %d" % (geno1, geno2, cpos))
    max1 = self.get_max_prob(gen_vals1, probidx)
    max2 = self.get_max_prob(gen_vals2, probidx)
    if max2 > max1:
      #print "G2: %s vs %s" % (geno1, geno2)
      return geno2 
    #print "G1: %s vs %s" % (geno1, geno2)
    return geno1

  def call_geno_for_threshold(self, geno, probidx, threshold):
    gen_vals = geno.split(":")
    max_prob, max_idx =  self.get_max_prob(gen_vals, probidx)
    if max_prob < threshold:
      gen_vals[0] = "./."
    else:
      gen_vals[0] = self.geno_strings[max_idx]

    return ":".join(gen_vals)

  def get_max_prob(self, gen_vals, probidx=1):
    probs = gen_vals[probidx].split(",")
    max_prob = 0.0
    max_prob_idx = -9
    for i, prob in enumerate(probs):
      if float(prob) > max_prob:
        max_prob = float(prob)
        max_prob_idx = i
    return max_prob, max_prob_idx

  def check_concordancies(self, data_list, assays, chipval):
    hwe_values = [0.0] * len(data_list)
    maf_values = [0.0] * len(data_list)
    obs = [0.0] * 3 
    exp = [0.0] * 3 
    allele_ref_1 = ""
    allele_alt_1 = ""
    allele_ref_2 = ""
    allele_alt_2 = ""
    #print "CHECK_CONC:", len(data_list)
    for i, vcf_record in enumerate(data_list): 
      if len(vcf_record) > 0:
        vcfr = VCFrecord(vcf_record)
        probidx = vcfr.get_probidx()
        homref_count, het_count, homalt_count, nc_count, miss_count = self.vcfr.get_allele_counts_from_array(data)
        allele_a, allele_b = vcfr.get_alleles()
        if allele_ref_1 == "":
          allele_ref_1 = allele_a
          allele_alt_1 = allele_b
          # Add 1 to prevent 0-divide
          obs[0] = homref_count + 1 
          obs[1] = het_count + 1 
          obs[2] = homalt_count + 1 
        else:
          allele_ref_2 = allele_a
          allele_alt_2 = allele_b
          exp[0] = homref_count + 1 
          exp[1] = het_count + 1 
          exp[2] = homalt_count + 1 
          if (allele_ref_1 != allele_ref_2) or (allele_alt_1 != allele_alt_2):
            varid = vcfr.get_varid(data)
            posn = vcfr.get_posn(data)
            self.allele_discord_count += 1
            logging.info("Allele discordancy: assay1=%s, assay2=%s, varid=%s, posn=%d, ref1=%s, alt1=%s, ref2=%s, alt2=%s", assays[0], assays[i], varid, int(posn), allele_ref_1, allele_alt_1, allele_ref_2, allele_alt_2)
            #print "Allele discord"
            return False

          chi_stat, chi_p_value = chisquare(obs, f_exp=exp)
          varid = vcfr.get_varid(data)
          posn = vcfr.get_posn(data)
          if chi_p_value < chipval:
            self.chisq_count += 1
            logging.info("CHI SQ test REJECT: assay1=%s, assay2=%s, varid=%s, posn=%d, chistat=%f, p_val=%e, chipval=%e, obs=%s, exp=%s, at %d", assays[0], assays[i], varid, int(posn), chi_stat, chi_p_value, chipval, str(obs), str(exp), i)
            #print "CHISQ discord"
            return False
          logging.info("CHI SQ test OK: assay1=%s, assay2=%s, varid=%s, posn=%d, chistat=%f, p_val=%e, obs=%s, exp=%s, at %d", assays[0], assays[i], varid, int(posn), chi_stat, chi_p_value, str(obs), str(exp), i)

    return True

  def check_maf_concordancies(self, data_list, maf_delta=0.3):
    maf_values = [0.0] * len(data_list)
