import os
import pymongo
import logging
from flask import Flask
from config import basedir
from config import LOGFILE
from models import DataStore

app = Flask(__name__)
app.config.from_object('config')

FORMAT = '%(asctime)-15s %(message)s'
logging.basicConfig(format=FORMAT, filename=LOGFILE ,level=logging.DEBUG)

logging.info("START")
#try:
ds = DataStore()
#except:
#  print "Unexpected error:", sys.exc_info()[0]

from app import views, models
