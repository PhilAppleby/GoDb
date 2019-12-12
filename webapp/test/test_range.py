import pymongo
import time
from zipfile import ZipFile
from zipfile import ZIP_DEFLATED
from optparse import OptionParser

from models_test import DataStore


# execution flow starts here
#
parser = OptionParser()
parser.add_option("-c", "--chr", dest="chr",
    help="chromosome", metavar="STR")
parser.add_option("-s", "--start", dest="start",
    help="start position", metavar="STR")
parser.add_option("-e", "--end", dest="end",
    help="end position", metavar="STR")
(options, args) = parser.parse_args()

connection = pymongo.MongoClient("mongodb://localhost")
#db = connection.gwas_metadata_2
db = connection.genomicsdb_test

ds = DataStore(db, 'akh')
start_time = time.time()

print "Extent:", int(options.end) - int(options.start)

print "range - get range data"
(sample_return_data, snp_return_data, msg) = ds.get_range_data(options.chr, options.start, options.end, 0.9, "combine", "csv")
print "Make zip data (r)"
zipfilename = "test_range_results.zip"
body = ds.make_zipfile(sample_return_data, snp_return_data, '.', zipfilename)
print "END:", time.time() - start_time, "seconds"

