#!/bin/sh
source /Users/phila/devel/godb/cfg/godb.cfg
ASSAYTYPE=$1
CHR=$2

mv ${DBDATADIR}/${ASSAYTYPE}/chr${CHR}.vcf.gz ${DBDATADIR}/${ASSAYTYPE}/chr${CHR}.vcf_old.gz

gzcat ${DBDATADIR}/${ASSAYTYPE}/chr${CHR}.vcf_old.gz | python ${PYUDIR}/convert_vcf_samples.py \
	--convfile=${DBDATADIR}/samples/prochi_convert0.csv \
	> ${DBDATADIR}/${ASSAYTYPE}/chr${CHR}.vcf

bgzip ${DBDATADIR}/${ASSAYTYPE}/chr${CHR}.vcf

tabix -pvcf -f ${DBDATADIR}/${ASSAYTYPE}/chr${CHR}.vcf.gz
