#!/bin/sh
export CONFFILE=$1
source ${CONFFILE}
date
export CHR=$2
${ZCATEXEC} ${DBDATADIR}/chr${CHR}.vcf.gz | python ${PYLDIR}/load_variants_from_vcf.py --assaytype=${ASSAYTYPE}
