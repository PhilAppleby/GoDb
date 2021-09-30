#!/bin/sh
export CONFFILE=$1
source ${CONFFILE}
python ${PYUDIR}/correct_and_flip.py \
	--csvfile=${EXTDATADIR}/enrique/ParaThyroidSNPS_GDctrlrsidupdated.csv \
	--snpsummary=${EXTDATADIR}/enrique/parathyoridGDCTRLfreq.afreq.csv \
	--dbsnpfile=${DBSNPDIR}/allchr_numeric_annotated.vcf.gz \
	--logfile=${LOGDIR}/correct_and_flip_gdtdtrls.log \
	> ${EXTDATADIR}/enrique/ParaThyroidSNPS_GDctrlrsidupdated_flipped.csv
