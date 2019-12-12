#!/bin/sh
export CONFFILE=$1
source ${CONFFILE}
perl ${PLLDIR}/process_geno_file_dir.pl --tmplt="%s/load_variants_from_single_info_file.sh %s %s" \
	--fsfx="info" --bindir=${SHLDIR} --dirname=${DBDATADIR} \
	--cfgfile=${CONFFILE} --printonly=$2 
