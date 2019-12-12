#!/bin/sh
export CONFFILE=$1
source ${CONFFILE}
perl ${PLLDIR}/process_geno_file_dir.pl --tmplt="%s/load_variants_from_single_vcf_file.sh %s %s" \
	--fsfx="vcf.gz" --bindir=${SHLDIR} --dirname=${DBDATADIR} \
	--cfgfile=${CONFFILE} --printonly=$2 
