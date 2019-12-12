import logging
# Helper methods for parsing VCF records - very basic just now
class VCFrecord():
  def __init__(self):
    self.first_genotype_idx = 9
    self.chr_idx = 0
    self.posn_idx = 1
    self.rsid_idx = 2
    self.ref_idx = 3
    self.alt_idx = 4
    self.qual_idx = 5
    self.filter_idx = 5
    self.info_idx = 7
    self.fmt_idx = 8
    self.calls = ["0/0", "0/1", "1/1", "./."]
    self.icalls = [0, 1, 2, -9]

  def get_first_genotype_index(self):
    """
    Return index of stat of genotype array
    """
    return self.first_genotype_idx

  def get_data_array(self, full_record):
    """
    The argument is an unsplit VCF record
    """
    return (full_record.split('\t'))

  def get_chr(self, full_record):
    """
    The argument is an unsplit VCF record
    """
    return (full_record.split('\t')[self.chr_idx])

  def get_chr_from_array(self, data):
    """
    The argument is a split VCF record
    """
    return (data[self.chr_idx])

  def get_varid(self, full_record):
    """
    The argument is an unsplit VCF record
    """
    return (full_record.split('\t')[self.rsid_idx])

  def get_varid_from_array(self, data):
    """
    The argument is a split VCF record
    """
    return (data[self.rsid_idx])
  
  def get_varid_ukb(self, full_record):
    """
    The argument is an unsplit VCF record
    """
    varid = full_record.split('\t')[self.rsid_idx]
    return (varid.split(',')[0])

  def get_varid_ukb_from_array(self, data):
    """
    The argument is a split VCF record
    """
    varid = data[self.rsid_idx]
    return (varid.split(',')[0])

  def get_alleles(self, full_record):
    """
    The argument is an unsplit VCF record
    """
    data = full_record.split('\t')
    return (data[self.ref_idx], data[self.alt_idx])

  def get_alleles_from_array(self, data):
    """
    The argument is a split VCF record
    """
    return (data[self.ref_idx], data[self.alt_idx])

  def get_posn(self, full_record):
    """
    The argument is an unsplit VCF record
    """
    return (full_record.split('\t')[self.posn_idx])

  def get_posn_from_array(self, data):
    """
    The argument is a split VCF record
    """
    return (data[self.posn_idx])

  def get_posn_as_int(self, full_record):
    """
    The argument is an unsplit VCF record
    """
    return (int(full_record.split('\t')[self.posn_idx]))

  def get_posn_from_array_as_int(self, data):
    """
    The argument is a split VCF record
    """
    return (int(data[self.posn_idx]))

  def get_info(self, full_record):
    """
    The argument is an unsplit VCF record
    """
    data = full_record.split('\t')
    return (data[self.info_idx].split(';'))

  def get_info_from_array(self, data):
    """
    The argument is a split VCF record
    """
    return (data[self.info_idx].split(';'))

  def set_info_in_array(self, data, info):
    """
    The argument is a split VCF record
    """
    data[self.info_idx] = info
    return data

  def set_varid_in_array(self, data, varid):
    """
    The arguments is a split VCF record and a variant id to be inserted
    """
    data[self.rsid_idx] = varid
    return data

  def parse_info(self, info_data):
    info = {}
    for info_elem in info_data:
      if "=" in info_elem:
        key_val = info_elem.split("=")
        info[key_val[0]]=key_val[1]
    return info

  def get_info_value_from_array(self, data, key):
    """
    The argument is a split VCF record
    """
    info = self.parse_info(self.get_info_from_array(data))
    #logging.info("info data: %s", str(info))
    if key in info:
      return info[key]
    else:
      return None

  def get_fmts(self, full_record):
    """
    The argument is an unsplit VCF record
    """
    data = full_record.split('\t')
    return (data[self.fmt_idx].split(':'))

  def get_fmts_from_array(self, data):
    """
    The argument is a split VCF record
    """
    #logging.info("get_fmts_from_array: %s", str(data[self.fmt_idx]))
    return (data[self.fmt_idx].split(':'))

  def get_fmt_indices(self, fmts, req_fmts):
    """
    Arguments are 2 lists of strings
    """
    idxs = {}
    for idx, fmt in fmts:
      if fmt in req_fmts:
        idxs[fmt] = idx
    return idxs

  def get_fmt_idx(self, fmts, req_fmt):
    """
    Arguments are 2 lists of strings
    """
    for i, fmt in enumerate(fmts):
      if fmt == req_fmt:
        return i
    return -9

  def get_probidx_from_array(self, data):
    return(self.get_fmt_idx(self.get_fmts_from_array(data), "GP"))

  def set_fmts_in_array(self, data, fmts):
    """
    The argument are a split VCF record and an array of formats
    """
    data[self.fmt_idx] = ":".join(fmts)
    return data

  def get_prfx_sfx(self, full_record):
    """
    The argument is an unsplit VCF record
    """
    rec = full_record.split('\t')
    return (rec[:self.first_genotype_idx], rec[self.first_genotype_idx:])

  def get_prfx_sfx_from_array(self, data):
    """
    The argument is an split VCF record
    """
    return (data[:self.first_genotype_idx], data[self.first_genotype_idx:])

  def get_call_rate(self, full_record):
    """
    """
    prfx,sfx = self.get_prfx_sfx(full_record)
    num_genotypes = len(sfx)
    call_count = 0
    for geno in sfx:
      try:
        geno_data = geno.split(":")
      except:
        continue
      if geno_data[0] != "./.":
        call_count += 1

    return (float(call_count) / num_genotypes)

  def get_call_rate_from_array(self, data):
    """
    """
    prfx,sfx = self.get_prfx_sfx_from_array(data)
    num_genotypes = len(sfx)
    call_count = 0
    for geno in sfx:
      try:
        geno_data = geno.split(":")
      except:
        continue
      if geno_data[0] != "./.":
        call_count += 1

    return (float(call_count) / num_genotypes)

  def get_genotype_counts(self, full_record):
    """
    """
    geno_counts = {}
    called_samp_count = 0
    samp_count = 0

    prfx, sfx = self.get_prfx_sfx(full_record)

    for geno in sfx:
      geno_data = geno.split(":")
      geno_counts[geno_data[0]] = geno_counts.get(varcall, 0) + 1
      samp_count += 1
      if geno_data[0] != "./.":
        called_samp_count += 1

    return geno_counts, called_samp_count, samp_count

  def get_genotype_counts_from_array(self, data):
    """
    """
    geno_counts = {}
    called_samp_count = 0
    samp_count = 0

    prfx, sfx = self.get_prfx_sfx_from_array(data)

    for geno in sfx:
      geno_data = geno.split(":")
      geno_counts[geno_data[0]] = geno_counts.get(geno_data[0], 0) + 1
      samp_count += 1
      if geno_data[0] != "./.":
        called_samp_count += 1

    return geno_counts, called_samp_count, samp_count

  def get_allele_counts_from_array(self, data, threshold=0.9):
    """
    """
    homref_count = 0
    het_count = 0
    homalt_count = 0
    nc_count = 0
    miss_count = 0 # can be "." or ""

    prfx, sfx = self.get_prfx_sfx_from_array(data)
    probidx = self.get_probidx_from_array(data)
    if probidx == -9:
      probidx = 1

    for geno in sfx:
      if geno == "." or geno == "":
        miss_count += 1
      else:
        geno_data = geno.split(":")
        probVals = geno_data[probidx]
        (call, ccode, maxprob) = self.get_call(probVals.split(","), threshold)
        if call == "0/0":
          homref_count += 1
        elif call == "1/1":
          homalt_count += 1
        elif call == "./.":
          nc_count += 1
        else:
          het_count += 1

    return homref_count, het_count, homalt_count, nc_count, miss_count
    
  def get_call(self, probs, threshold):
    max_prob = 0.0
    max_idx = 3

    #print "Probs:", probs
    for idx, prob in enumerate(probs):
      if float(prob) > max_prob:
        max_prob = float(prob)
        max_idx = idx
 
    if (threshold !=0.0):
      #print 'threshold', threshold, max_prob
      if max_prob < threshold:
        #print 'LT threshold', threshold, max_prob
        max_idx = 3
 
    return (self.calls[max_idx], self.icalls[max_idx], max_prob)

  def get_call_from_probs(self, gen_vals, probidx=1, threshold=0.9):
    probs = gen_vals[probidx].split(",")
    max_prob = 0.0 
    max_idx = 3 
    
    for idx, prob in enumerate(probs):
      if float(prob) > max_prob:
        max_prob = float(prob)
        max_idx = idx 

    if (float(threshold) !=0.0):
      if max_prob < float(threshold):
        max_idx = 3 
    
    return (self.calls[max_idx], self.icalls[max_idx], max_prob, max_idx)


