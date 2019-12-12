import pymongo
import time
from zipfile import ZipFile
from zipfile import ZIP_DEFLATED
from optparse import OptionParser

from models import DataStore


# execution flow starts here
#
parser = OptionParser()
parser.add_option("-r", "--rsid", dest="rsid",
    help="dbSNP rsid", metavar="STR")
parser.add_option("-p", "--prochi", dest="prochi",
    help="prochi", metavar="STR")
(options, args) = parser.parse_args()

connection = pymongo.MongoClient("mongodb://localhost")
db = connection.gwas_metadata_2
#db = connection.gwas_metadata

ds = DataStore(db, 'akh')
start_time = time.time()

ds.get_rsid_prochi_data(options.rsid, options.prochi, 0.9, "csv_calls")
connection.close()

print "END:", time.time() - start_time, "seconds"
