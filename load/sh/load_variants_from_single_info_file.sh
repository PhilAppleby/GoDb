#!/bin/sh
export CONFFILE=$1
source ${CONFFILE}
date
export CHR=$2
cat ${DBDATADIR}/chr${CHR}.info | python ${PYLDIR}/load_variants_from_info.py --chr=${CHR} --assaytype=${ASSAYTYPE}
