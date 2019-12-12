#!/bin/sh
export CONFFILE=$1
source ${CONFFILE}
python ${PYLDIR}/load_filepaths.py --prfx=${DBDATAPRFX} --sfx=${DBDATASFX} --assaytype=${ASSAYTYPE}
