#!/bin/sh
export PYLDIR=/Users/phila/devel/godb/load/py
export DATADIR=/Users/phila/data/
cat ${DATADIR}/snp_results/rs7412_local/$1_samples.csv | python ${PYLDIR}/compare_csv_genotypes.py \
	--convfile=${DATADIR}/godb/samples/all_godarts_prochis_convert.csv \
       	--genofile=${DATADIR}/snp_results/rs7412_remote/$1_samples.csv
