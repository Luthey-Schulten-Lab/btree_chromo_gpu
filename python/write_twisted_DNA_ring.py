# =============================================
# Author: Benjamin R. Gilbert
# Email: brg4@illinois.edu
# =============================================

import numpy as np

# import time
# import os

# import importlib
# import sys

# sys.path.insert(1,'./system_preparation')
from system_preparation import DNA_rings as D_r
# imp.reload(D_r)

out_dir = '/home/andrew/Data/btree_chromo/twist_testing/'
out_label = 'initial_ring'

# r = 17.0
# R = 600.0
# c = 0.91
r = 17.0
R = 275.0
c = 1.0
ntwists = 0
n, x, q = D_r.initialize_ring(r,R,c,ntwists,1.0,1.0)
# r = 17.0
# R = 300.0
# c = 14.0
#n, x, q = D_r.initialize_toroid(r,R,c,0,1.0,1.0,100.0,40)
# r = 17.0
# R = 150.0
# c = 3.0
# n, x, q = D_r.initialize_trefoil(r,R,c,0)

print(x.shape)

chromo_file = out_dir + 'chromo_' + out_label + '.txt'

with open(chromo_file,'w') as f:
    f.write('size='+str(n))

coord_file = out_dir + 'x_mono_' + out_label + '.bin'

with open(coord_file,'wb') as f:
    x.tofile(f)

quat_file = out_dir + 'q_mono_' + out_label + '.bin'

with open(quat_file,'wb') as f:
    q.tofile(f)

with open(quat_file,'rb') as f:
    newq = np.fromfile(quat_file, dtype=np.double)
    newq = np.reshape(newq, q.shape, 'C')
    print(newq)