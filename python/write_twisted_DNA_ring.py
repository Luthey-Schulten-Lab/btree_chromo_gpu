# =============================================
# Author: Benjamin R. Gilbert
# Email: brg4@illinois.edu
# =============================================

import numpy as np

import time
import os

import imp
import sys

sys.path.insert(1,'./system_preparation')
import DNA_rings as D_r
imp.reload(D_r)

out_dir = '/home/ben/Data/btree_chromo/persistence_length_testing/'
out_label = 'test'

r = 17.0
R = 1000.0
c = 1.8

n, x, q = D_r.initialize_ring(r,R,c,20,1.0,1.0)

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

