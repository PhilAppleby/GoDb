#!/bin/sh
export CONFFILE=$1
source ${CONFFILE}
python ${PYUDIR}/correct_and_flip.py \
	--csvfile=${EXTDATADIR}/enrique/ParaThyroidSNPS_GoSHAREfreeze1.csv \
	--snpsummary=${EXTDATADIR}/enrique/parathyoridmergedgoshare.afreq.csv \
	--dbsnpfile=${DBSNPDIR}/allchr_numeric_annotated.vcf.gz \
	--logfile=${LOGDIR}/correct_and_flip_goshare.log \
	> ${EXTDATADIR}/enrique/ParaThyroidSNPS_GoSHAREfreeze1_flipped.csv