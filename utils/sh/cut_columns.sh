#!/bin/sh
export CONFFILE=$1
source ${CONFFILE}
python ${PYUDIR}/cut_csv_geno_file.py \
	--colnames=$2 \
	--logfile=${LOGDIR}/cut_csv_geno.log \
	--csvfile=${EXTDATADIR}/enrique/$3 \
	> ${EXTDATADIR}/enrique/$4
