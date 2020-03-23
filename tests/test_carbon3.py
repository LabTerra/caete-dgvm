
import sys
sys.path.insert(0, '../src/')
from caete_module import global_par as gp
from caete_module import soil_dec as carbon3
from caete_module import soil_pools as funcs
import csv
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def zeroes(*args):
    return np.zeros(shape=args, dtype=np.float32)


# start inputs
lnr = zeroes(6) + 0.0001


#    subroutine carbon3(tsoil, water_sat, leaf_litter, coarse_wd,&
#                     &        root_litter, lnr, cl, cs, &
#                     &         snr_in,
#
# avail_p,
# inorg_n,
# inorg_p,&
# sorbed_p
# cl_out,
# cs_out,
# snr,
# hr
header = ['hr', 'avail-p', 'inorg_n', 'inorgp', 'sorbed_p', 'cl1',
          'cl2', 'cs1', 'cs2', 'snr1-cl1_n2c', 'snr2-cl2_n2c', 'snr3', 'snr4', 'snr5', 'snr6', 'snr7', 'snr8']


with open('carbon3_test.csv', 'w') as fh:
    CSV_WRITER = csv.writer(fh, delimiter=',')
    CSV_WRITER.writerow(header)
    for x in range(300000):
        # UPDATE_ input VARIABLES
        snr_in = funcs.sp_snr
        avail_p = funcs.sp_available_p
        inorg_n = funcs.sp_available_n
        inorg_p = funcs.sp_in_p
        sorbed_p = funcs.sp_so_p
        cl = funcs.sp_csoil[0:2]
        cs = funcs.sp_csoil[2:]
        df = carbon3.carbon3(23, 0.5, 1, 1, 1, lnr, cl, cs,
                             snr_in, avail_p, inorg_n, inorg_p, sorbed_p)
        # print(df)
        line = [df[-1],
                df[0],
                df[1],
                df[2],
                df[3],
                df[4][0],
                df[4][1],
                df[5][0],
                df[5][1],
                df[6][0],
                df[6][1],
                df[6][2],
                df[6][3],
                df[6][4],
                df[6][5],
                df[6][6],
                df[6][7]]
        CSV_WRITER.writerow(line)
        # update globl pools

        funcs.sp_snr = df[6][:]
        funcs.sp_csoil[0:2] = df[4][:]
        funcs.sp_csoil[2:] = df[5][:]
        funcs.sp_available_p = df[0]
        funcs.sp_available_n = df[1]
        funcs.sp_in_p = df[2]
        funcs.sp_so_p = df[3]

data = pd.read_csv('carbon3_test.csv')
os.system('rm -rf carbon3_test.csv')
