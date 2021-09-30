#!/bin/sh
./cut_columns_file_pairs.sh ../../cfg/godb.cfg $1 \
	ParaThyroidSNPS_GoSHAREfreeze1_flipped.csv ParaThyroidSNPS_GoSHAREfreeze1_columns.csv \
	cut_goshare_file.csv cut_goshare_file_columns.csv
./cut_columns_file_pairs.sh ../../cfg/godb.cfg $1 \
	ParaThyroidSNPS_GDctrlrsidupdated_flipped.csv ParaThyroidSNPS_GDctrlrsidupdated_columns.csv \
	cut_godartsctrls_file.csv cut_godartsctrl_columns.csv
