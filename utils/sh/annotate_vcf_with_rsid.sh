#!/bin/sh
export CONFFILE=$1
source ${CONFFILE}
CHR=`printf "%.2d" ${2}`
gzcat ${MIDATADIR}/${CHR}_dose2020_rsid.vcf.gz | python ${PYLDIR}/annotate_vcf_with_rsid.py --chrom=${CHR} \
	--chrommap=${DBSNPDIR}/chrom_map.txt \
	--logfile=${LOGDIR}/${CHR}_annot.log \
	--dbsnpfile=${DBSNPDIR}/GCF_000001405.25_Chr${CHR}_rsids.vcf.gz > \
	${DBSNPDIR}/chr${CHR}_annotated.vcf
