import pymongo
import time
from optparse import OptionParser

from models import DataStore

start_time = time.time()

parser = OptionParser()
parser.add_option("-s", "--rsid", dest="rsid",
  help="rsid, eg rs7412", metavar="STRING")
(options, args) = parser.parse_args()

connection = pymongo.MongoClient("mongodb://localhost")
db = connection.metadata

ds = DataStore(db, 'affy')

print "STUDY NAME", ds.get_studyname()
marker_data = ds.get_marker_summary_probs(options.rsid)
print marker_data["TOTALS"]
connection.close()

print "END:", time.time() - start_time, "seconds"

