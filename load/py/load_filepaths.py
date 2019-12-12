#
# Add a record to the mongodb genomicsdb_test filepaths collection
# 
import time
import os, sys
from optparse import OptionParser
from godb import GoDb

start_time = time.time()
filedocs = []

def main(options):
  try:
    godb = GoDb()
  except:
    print "Unexpected error:", sys.exc_info()[0]
    exit()

  filelist = []
  dirname = "/" + options.prfx + "/" + options.sfx + "/"

  for filename in os.listdir(dirname):
    if filename.endswith('.vcf.gz'):
      filelist.append(filename)

  godb.add_filepath_detail(options.assaytype, options.prfx, options.sfx, filelist)
#
# execution flow starts here
#
parser = OptionParser()
parser.add_option("-s", "--sfx", dest="sfx",
   help="Directory suffix for .vcf.gz files", metavar="DIR")
parser.add_option("-p", "--prfx", dest="prfx",
   help="Directory prefix for .vcf.gz files", metavar="DIR")
parser.add_option("-a", "--assaytype", dest="assaytype",
   help="Assay type (illumina, affy, etc)", metavar="STR")

start_time = time.time()

(options, args) = parser.parse_args()

rec_count = main(options)
print "END:", time.time() - start_time, "seconds", rec_count

