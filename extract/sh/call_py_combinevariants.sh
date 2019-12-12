#!/bin/sh
# Auto build gene scores
# parameter $1 is subpath (can be project or project/subproject)
# Files are sources as follows:
# .../combined_samples.csv - downloaded from the GoDARTS_GDb page
# .../standardised_weights.csv - snp data and weights from other research, in a standard format.
# 
source ${HOME}/devel/godb/cfg/godb.cfg
#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------
python ${PYEDIR}/combine_variants.py --check=N --snpfile=$1  --prfx=${DBDATAPRFX} --logfile=${LOGDIR}/merge_extract.log > ${COMBODATADIR}/genotypes.txt
grep -v OVERLAP ${COMBODATADIR}/genotypes.txt > ${COMBODATADIR}/genotypes.vcf
