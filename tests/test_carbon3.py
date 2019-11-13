import csv
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import sys

sys.path.insert(0, '../src/')

from caete_module import soil_dec as funcs
from caete_module import global_par as gp

def zeroes(*args):
    return np.zeros(shape=args, dtype=np.float32)


cl = zeroes(2) + 0.01
cs = zeroes(3) + 0.01
nmin = zeroes(1) + 0.02
plab = zeroes(1) + 0.003
lnr = zeroes(6) + np.random.random(6,) * 0.001
header = ['hr','cl1','cl2','cs1','cs2','cs3']

#   subroutine carbon3(tsoil,leaf_l,cwd,root_l,lnr,cl,cs,cl_out,cs_out,hr)
# tsoil, leaf_l, cwd, root_l, lnr, cl, cs, cl_out, cs_out, hr
with open('carbon3_test.csv', 'w') as fh:
    CSV_WRITER = csv.writer(fh, delimiter=',')
    CSV_WRITER.writerow(header)
    for x in range(100000):
        df = funcs.carbon3(20, 0.1, 0.1, 0.1, lnr, cl, cs)
        line = [df[-1],
                df[0][0],
                df[0][1],
                df[1][0],
                df[1][1],
                df[1][2]]
        CSV_WRITER.writerow(line)
        cl = df[0][:]
        cs = df[1][:]


data = pd.read_csv('carbon3_test.csv')
os.system('rm -rf carbon3_test.csv')
