#!/bin/sh
# Note that in a production system this should be run as a gateway service
# Example start / stop commands
# "sudo systemctl start genomicsdb"
# "sudo systemctl stop genomicsdb"
date
CONFFILE=$1
source $CONFFILE
#python ${PYWDIR}/run.py > ${LOGDIR}/godb_server.log 2>&1 &
python ${PYWDIR}/run.py 

