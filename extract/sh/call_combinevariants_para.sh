#!/bin/sh
#---------------------------------------------------------------------------------------------------------------------
# Data extract and produce combined genotype records the GoDb database
# parameter $1 is full path to a SNP list
#---------------------------------------------------------------------------------------------------------------------
# 
source ${HOME}/devel/godb/cfg/godb.cfg
#---------------------------------------------------------------------------------------------------------------------
# Get data from the GoDb 
#---------------------------------------------------------------------------------------------------------------------
#echo ${ASSAYTYPES}
${GOBIN}/vcombine -rsfile=$1 \
	-vcfprfx=${DBDATAPRFX} -threshold=0.9 \
	-assaytypes=${ASSAYTYPES} \
	> ${COMBODATADIR}/merge_output.txt
#---------------------------------------------------------------------------------------------------------------------
# Separate output files
#---------------------------------------------------------------------------------------------------------------------
#grep -w combined ${COMBODATADIR}/merge_output.txt | grep -v METRICS | cut -f2- > ${COMBODATADIR}/combined.vcf
