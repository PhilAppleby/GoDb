# -*- coding: utf-8 -*-
import os
basedir = os.environ['DATADIR']

CSRF_ENABLED = True
SECRET_KEY = 'you-will-never-guess'

EXTERNAL_CONNECT_PORT = 80
#UPLOAD_DIR = basedir + '/' + 'uploads'
UPLOAD_DIR = '/tmp'
LOGFILE = basedir + '/logs/' + 'gdb_search.log'
ALLOWED_EXTENSIONS = set(['txt'])
RSLISTMAX = 200
