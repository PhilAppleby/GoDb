# Normalise prescribing data 
# 
import time
import os, sys

start_time = time.time()

def main():

  hdrline = sys.stdin.readline()
  hdr = hdrline.strip().split(",")
  for i, colname in enumerate(hdr):
    #print i, colname
    print colname

# execution flow starts here
#
main()

