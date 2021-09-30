class GRShelper():
  def __init__(self):
    self.opp = {
      'G':'C',
      'A':'T',
      'C':'G',
      'T':'A',
    }

  def load_snp_weight_data_common(self, fh):
    """
    Load SNP weight data
    (Often derived from other studies)
    EAF (effect allele frequency) and RAF (risk allele frequency)
    are used interchangeably
    """
    rsids = []
    snp_wts = {}
    eff_alleles = {}
    rafs = {}
    hdr = fh.readline().strip().split(',')
    rsid_col = self.get_col_num(hdr, 'rsid')
    ea_col = self.get_col_num(hdr, 'ea')
    eaf_col = self.get_col_num(hdr, 'eaf')
    wgt_col = self.get_col_num(hdr, 'wgt')

    for line in fh:
      data=line.strip().split(',')
      rsids.append(data[rsid_col])
      eff_alleles[data[rsid_col]] = data[ea_col].upper()
      eaf = data[eaf_col].upper()
      if eaf == "":
        eaf = 0.0
      rafs[data[rsid_col]] = eaf
      wgt = data[wgt_col]
      if wgt == "":
        wgt = 0.0
      snp_wts[data[rsid_col]] = wgt

    return rsids, snp_wts, eff_alleles, rafs

  def get_col_num(self, hdr_arr, col_txt):
    for i, elem in enumerate(hdr_arr):
      if col_txt == elem:
        return i

    raise ValueError("Col %s not found" % (col_txt))

  def load_snp_summary_data(self, fh):
    """
    SNP summary data is derived from GoDARTS_GDb data
    """
    alleles = {}

    for line in fh:
      data=line.strip().split(',')
      alleles[data[0]] = data[1] + "/" + data[2]

    return alleles

  def isCompatible(self, wtRiskAllele, wtRaf, locRiskAllele, locNonRiskAllele, locRaf):
    """
    Determine if we should compare local data with incoming snp weight data
    """
    if wtRiskAllele == locRiskAllele:
      return True

    if wtRiskAllele == locNonRiskAllele:
      return True

    if wtRiskAllele == self.opp[locRiskAllele]:
      return True

    if wtRiskAllele == self.opp[locNonRiskAllele]:
      return True

    # 0.48 and 0.52 are somewhat arbitrary ...
    if (wtRaf > 0.48 and wtRaf < 0.52):
      return False

    if (locRaf > 0.48 and locRaf < 0.52):
      return False

    return False

  def notEqual(self, wtRa, locRa):
    """
    Compare alleles
    """
    if (wtRa == locRa) or (wtRa == self.opp[locRa]):
      return False
#   inequality
    return True

  def should_flip(self, wtRaf, locRaf):
    """
    """
    if wtRaf < 0.5 and locRaf < 0.5:
      return False
    if wtRaf > 0.5 and locRaf > 0.5:
      return False

    return True

  def should_flip_full(self, rsid, wtRa, wtRaf, locNra, locRa, locRaf):
    """
    Check allele letter values and rafs where there are potential
    (AT, CG) strand conflicts in which case a comparison of Allelle frequencies
    is used - if they're different side of 0.5 don't flip
    """
    if wtRa == locRa:
      print "NOFLIP,SAME,%s,%s,%s" % (rsid, wtRa, locRa)
      return False

    if wtRa == self.opp[locRa]:
      if locNra == self.opp[locRa]:
        if not self.should_flip(wtRaf, locRaf):
          print "NOFLIP,STRAND_MAF,%s,%s,%s" % (rsid, wtRa, locRa)
          return False
      else:
        print "NOFLIP,OPP,%s,%s,%s" % (rsid, wtRa, locRa)
        return False

    return True
