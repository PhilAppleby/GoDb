#!/bin/sh
export CONFFILE=$1
source ${CONFFILE}
tabix -h ${DBDATADIR}/chr19.vcf.gz 0:0-0 | python ${PYLDIR}/load_samples.py --assaytype=${ASSAYTYPE}

