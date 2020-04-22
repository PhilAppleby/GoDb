#!/bin/sh
CONFFILE=$1
echo ${CONFFILE}
source ${CONFFILE}
date
cat ${GMAPDIR}/refFlat.txt | python ${PYLDIR}/load_gene_map.py
