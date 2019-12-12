from vcfrecord import VCFrecord

class Mafhelper():

  def get_maf_and_cr(self, vcfr):
    """
    Violates the do one thing rule, but means only 1 pass over the genotype_array
    """
# Get some metrics from the genotype_array
    ref, alt = vcfr.get_alleles()
    homref_count, het_count, homalt_count, nc_count, miss_count = vcfr.get_allele_counts()
    call_count = homref_count + het_count + homalt_count
    cr = float(call_count)/(call_count + nc_count)
# MaF
    maf, ma = self.maf(het_count, homref_count, ref, homalt_count, alt, nc_count)

    return maf, ma, cr

  def maf(self, obs_hets, obs_hom1, allele1, obs_hom2, allele2, num_samples):
    if obs_hom1 < 0 or obs_hom2 < 0 or obs_hets < 0:
      raise Exception("FATAL ERROR - MAF: Current genotype configuration (%s  %s %s) includes negative count" % (obs_hets, obs_hom1, obs_hom2))

    obs_hom_maj = obs_hom2 if obs_hom1 < obs_hom2 else obs_hom1
    obs_hom_min = obs_hom1 if obs_hom1 < obs_hom2 else obs_hom2
    ma = allele1 if obs_hom1 < obs_hom2 else allele2

    gcount = obs_hets + obs_hom_maj + obs_hom_min
    #gcount = num_samples

    minor_count = (obs_hom_min * 2) + (obs_hets)
    #print "MAF Numbers = ", gcount, minor_count

    try:
      return ("%.6f" % (float(minor_count) / (gcount * 2.0))), ma
    except:
      return 0.0, ma


  def af(self, obs_hets, obs_hom_interest, obs_hom_other, num_samples):
    gcount = obs_hets + obs_hom_interest + obs_hom_other

    allele_count = (obs_hom_interest * 2)  + (obs_hets)
    #print "AF Numbers = ", gcount, allele_count

    try:
      return ("%.6f" % (float(allele_count) / (gcount * 2.0)))
    except:
      return 0.0
