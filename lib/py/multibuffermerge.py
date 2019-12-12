from multimerge import Multimerge
# Helper methods for merging multiple record buffers
class Multibuffermerge(Multimerge):
  def __init__(self, hdr_cols):
    Multimerge.__init__(self, len(hdr_cols))
    self.prt_comments = []
    for i, hdr_list in enumerate(hdr_cols):
      #print "MB COLS", i, len(hdr_list), hdr_list
      self.hdr_cols[i] = hdr_list
    self._load_col_data()

