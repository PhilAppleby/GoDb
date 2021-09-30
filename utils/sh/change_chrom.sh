#!/bin/sh
export CONFFILE=$1
source ${CONFFILE}
gzcat ${DBSNPDIR}/GCF_000001405.25.gz | python ${PYLDIR}/change_chrom.py --chrommap=${DBSNPDIR}/chrom_map.txt \
	--logfile=${LOGDIR}/${CHR}_change.log >\
	${DBSNPDIR}/allchr_numeric_annotated.vcf
