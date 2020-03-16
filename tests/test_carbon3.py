
import sys
sys.path.insert(0, '../src/')
from caete_module import global_par as gp
from caete_module import soil_dec as funcs
import csv
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def zeroes(*args):
    return np.zeros(shape=args, dtype=np.float32)


cl = zeroes(2) + 0.01
cs = zeroes(2) + 0.01
nupt = 0.02
pupt = 0.01
lnr = zeroes(6) + np.random.random(6,) * 0.001
snr_in = zeroes(8)

avail_p = 0.0
inorg_n = 0.0
inorg_p = 0.0
sorbed_p = 0.0
cl = zeroes(2)
cs = zeroes(2)
snr = 0.0
hr = 0.0
avail_p, inorg_n, inorg_p, sorbed_p
header = ['hr', 'cl1', 'cl2', 'cs1', 'cs2', 'snr1',
          'snr2', 'snr3', 'snr4', 'snr5', 'snr6', 'snr7', 'snr8']


#    subroutine carbon3(tsoil, water_sat, leaf_litter, coarse_wd,&
#                     &        root_litter, lnr, cl, cs, &
#                     &         snr_in, avail_p, inorg_n, inorg_p,&
#                     &          sorbed_p, cl_out, cs_out, snr, hr

with open('carbon3_test.csv', 'w') as fh:
    CSV_WRITER = csv.writer(fh, delimiter=',')
    CSV_WRITER.writerow(header)
    for x in range(1):
        df = funcs.carbon3(23, 0.5, 20, 20, 20, lnr, cl, cs,
                           snr_in, avail_p, inorg_n, inorg_p, sorbed_p)
        print(df)
        # line = [df[-1],
        # df[0][0],
        # df[0][1],
        # df[1][0],
        # df[1][1],
        # df[2][0],
        # df[2][1],
        # df[2][2],
        # df[2][3],
        # df[2][4],
        # df[2][5],
        # df[2][6],
        # df[2][7]]
        # CSV_WRITER.writerow(line)
        # cl = df[0][:]
        # cs = df[1][:]


# data = pd.read_csv('carbon3_test.csv')
# os.system('rm -rf carbon3_test.csv')
