# ----------------------------------------------------------------
# Command line args:
# ----------------------------------------------------------------
#
import time
import pymongo
import re
import os, sys
from optparse import OptionParser

from hgdb import Hgdb # db helper methods

start_time = time.time()
output_data_size = 100
fa_rec_size = 50
genome_str_len = output_data_size * fa_rec_size

def main():

  try:
    connection = pymongo.MongoClient("mongodb://localhost")
  except IOError as e:
    print "I/O error({0}): {1}".format(e.errno, e.strerror)
    exit()
  except TypeError as e:
    print "Missing arguments ", e
    exit()
  except:
    print "Unexpected error:", sys.exc_info()[0]
    sys.exit()

  db = connection.hg
  hgdb = Hgdb(db, genome_str_len)

  count = 0
  glength = 0
  offset = 1
  output_total = 0
  output_data = []
  chromosome = 0

  for line in sys.stdin:
    line = line.strip()
    if line.startswith(">"):
        if len(output_data) > 0:
          out_str = "".join(output_data)
          hgdb.add_genome_data(chromosome, offset, out_str, start_time)
          output_data = []
        offset = 1
        line = line[1:]
        print line[0:3], line[3:], offset
        if line[3] == 'X':
          chromosome = 23
        elif line[3] == 'Y':
          chromosome = 24
        else:
          chromosome = int(line[3:])
    else:
      if len(output_data) >= output_data_size:
        out_str = "".join(output_data)
        output_total += 1
        hgdb.add_genome_data(chromosome, offset, out_str, start_time)
        #print offset, out_str, len(out_str)
        offset += len(out_str)
        #print "OFFSET %s -> %d (output_data_size = %d)" % (chromosome, offset, output_data_size)
        output_data = []
      count+= 1
      glength += len(line)
      output_data.append(line)

  out_str = "".join(output_data)
  hgdb.add_genome_data(chromosome, offset, out_str, start_time)
  hgdb.flush_genome_buffer(start_time)
  output_total += 1
  #print offset, out_str, len(out_str)
  return count, glength, output_total


# execution flow starts here
#
parser = OptionParser()

rec_count, gen_length, out_total = main()
print("END: {0:.5f} seconds stats={1},{2},{3},{4}".format(time.time() - start_time, rec_count, gen_length, out_total, genome_str_len))
