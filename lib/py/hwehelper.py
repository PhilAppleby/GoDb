class Hwehelper():
# Derived from code originally written in C by J, Wiggington:
# Available at: http://csg.sph.umich.edu//abecasis/Exact/snp_hwe.c
# The following comments were preserved:
#/*
#// This code implements an exact SNP test of Hardy-Weinberg Equilibrium as described in
#// Wigginton, JE, Cutler, DJ, and Abecasis, GR (2005) A Note on Exact Tests of 
#// Hardy-Weinberg Equilibrium. American Journal of Human Genetics. 76: 000 - 000  
#//
#// Written by Jan Wigginton
#*/
  def HWE_exact(self, obs_hets, obs_hom1, obs_hom2, num_samples):
    if obs_hom1 < 0 or obs_hom2 < 0 or obs_hets < 0:
      raise Exception("FATAL ERROR - Hwehelper.HWE_exact: Supplied genotype counts (%d  %d %d) contain a -ve amount" % (obs_hets, obs_hom1, obs_hom2))

    obs_homc = obs_hom2 if obs_hom1 < obs_hom2 else obs_hom1
    obs_homr = obs_hom1 if obs_hom1 < obs_hom2 else obs_hom2

    #print "hom minor: ", obs_homr
    #print "hom major: ", obs_homc
    #print "het      : ", obs_hets

    rare_copies = 2 * obs_homr + obs_hets
    common_copies = 2 * obs_homc + obs_hets
    #genotypes = obs_hets + obs_homc + obs_homr
    genotypes = num_samples

    #print rare_copies, common_copies, genotypes
    het_probs = [0.0] * (rare_copies + 1)

    #print "het_probs len at start: ", len(het_probs)
    #print "het_probs at start: ", het_probs
    #start at midpoint
    mid = rare_copies * (2 * genotypes - rare_copies) / (2 * genotypes)
    #print "mid(1):", mid

    #check to ensure that midpoint and rare alleles have same parity
    #print (rare_copies & 1) ^ (mid & 1)
    if (rare_copies & 1) ^ (mid & 1):
      mid += 1

    #print "mid(2):", mid
    curr_hets = mid
    curr_homr = (rare_copies - mid) / 2
    curr_homc = genotypes - curr_hets - curr_homr

    het_probs[mid] = 1.0
    sum = float(het_probs[mid])

    for curr_hets in xrange(mid, 1, -2):
      het_probs[curr_hets - 2] = het_probs[curr_hets] * curr_hets * (curr_hets - 1.0) / (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0))

      sum += het_probs[curr_hets - 2];

      # 2 fewer heterozygotes for next iteration -> add one rare, one common homozygote
      curr_homr += 1
      curr_homc += 1
    #print het_probs
    curr_hets = mid
    curr_homr = (rare_copies - mid) / 2
    curr_homc = genotypes - curr_hets - curr_homr

    for curr_hets in xrange(mid, rare_copies - 1, 2):

      het_probs[curr_hets + 2] = het_probs[curr_hets] * 4.0 * curr_homr * curr_homc / ((curr_hets + 2.0) * (curr_hets + 1.0))

      sum += het_probs[curr_hets + 2]

      #add 2 heterozygotes for next iteration -> subtract one rare, one common homozygote
      curr_homr -= 1
      curr_homc -= 1

    
    #print het_probs
    for i in xrange(0, rare_copies + 1):
      het_probs[i] /= sum

    #print het_probs
    #alternate p-value calculation for p_hi/p_lo
    p_hi = float(het_probs[obs_hets])
    for i in xrange(obs_hets, rare_copies+1):
      p_hi += het_probs[i]

    p_lo = float(het_probs[obs_hets])
    for i in xrange(obs_hets-1, -1, -1):
      p_lo += het_probs[i]

    p_hi_lo = 2.0 * p_hi if p_hi < p_lo else 2.0 * p_lo

    p_hwe = 0.0
    #  p-value calculation for p_hwe
    for i in xrange(0, rare_copies + 1):
      if het_probs[i] > het_probs[obs_hets]:
        continue;
      p_hwe += het_probs[i]

    p_hwe = 1.0 if p_hwe > 1.0 else p_hwe

    return "%.8f" % (p_hwe)

  def hwe_check(self, hets, hom1, hom2):
    minor_obs = 0
    major_obs = 0
    if hom1 > hom2:
      minor_obs = hom2
      major_obs = hom1
    else:
      minor_obs = hom1
      major_obs = hom2

    num_gt = hets + hom1 + hom2
    minor_freq = ((2 * minor_obs) + hets) / (num_gt * 2.0) 
    major_freq = ((2 * major_obs) + hets) / (num_gt * 2.0) 
    print (major_freq * major_freq)
    print (2 * major_freq * minor_freq)
    print (minor_freq * minor_freq)

    return((major_freq * major_freq) + (2 * major_freq * minor_freq) + (minor_freq * minor_freq))

  def hwe_simple(self, hets, hom1, hom2):
    minor_obs = 0
    major_obs = 0
    if hom1 > hom2:
      minor_obs = hom2
      major_obs = hom1
    else:
      minor_obs = hom1
      major_obs = hom2

    num_gt = hets + hom1 + hom2
