# Load the mongodb genemap collection
# 
import time
import os, sys
from godb import GoDb

start_time = time.time()
flush_at = 10000

def main():
  try:
    godb = GoDb()
  except:
    print "Unexpected error:", sys.exc_info()[0]
    exit()

  hdr = []
  hdrlen = 0
  count = 0
  for line in sys.stdin:
    line = line.strip()
    if (line.startswith('#')):
      pass
    else:
      godb.process_genemap_detail(line)
      count += 1 
      if (godb.get_genemap_len() >= flush_at):
        godb.flush_genemap_buff()
        print ".", time.time() - start_time, "seconds", count

  godb.flush_genemap_buff()
  print ""
  return count
#
# execution flow starts here
#
rec_count = main()
print "END:", time.time() - start_time, "seconds", rec_count
