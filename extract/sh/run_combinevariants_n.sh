#!/bin/sh
#---------------------------------------------------------------------------------------------------------------------
# Data extract and produce combined genotype records the GoDb database
# parameter $1 is full path to a SNP list
#---------------------------------------------------------------------------------------------------------------------
# 
source ${HOME}/devel/godb/cfg/godb.cfg


for i in {1..100}
do
	${SHEDIR}/call_combinevariants.sh $1
done
