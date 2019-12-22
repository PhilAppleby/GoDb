#!/bin/sh
# Extract combined SNP records from the GoDb
# 
source ${HOME}/devel/godb/cfg/godb.cfg
#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------
python ${PYEDIR}/combine_variants.py --check=N --snpfile=$1  --prfx=${DBDATAPRFX} --logfile=${LOGDIR}/merge_extract.log > ${COMBODATADIR}/genotypes.txt
grep -v OVERLAP ${COMBODATADIR}/genotypes.txt > ${COMBODATADIR}/genotypes.vcf
