import sys
sys.path.append('C:/Users/maple/Desktop/代码/iDIRECT')
import idirect as idir
import file_handler as fh
import net_handler as nh
from time import time
t0 = time()
G,n = fh.read_file_weighted_edges("C:/Users/maple/Desktop/代码/iDIRECT/rmc.txt", t0)
S,err = idir.direct_association(G, t0=t0)
S2 = nh.merge(G, S)
St = fh.save_sorted_turple(S2, in_file="C:/Users/maple/Desktop/代码/iDIRECT/rmc.txt")
fh.save_file_weighted_edges(St, "C:/Users/maple/Desktop/代码/iDIRECT/rmc_res.txt")