# Load the mongodb metadata markers collection
#
import time
import os, sys

start_time = time.time()

def main():
  outlist = []
  for line in sys.stdin:
    outlist.append(line.strip())

  print ",".join(outlist)
  return
#
# execution flow starts here
#
main()
