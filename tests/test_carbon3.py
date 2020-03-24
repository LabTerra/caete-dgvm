
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


#    subroutine carbon3(tsoil, water_sat, leaf_litter, coarse_wd,&
#                     &        root_litter, lnr, cl, cs, &
#                     &         snr_in,
# avail_p
# inorg_n,
# inorg_p,
# sorbed_p
# cl_out,
# cs_out,
# snr,
# hr
# LNR [(l1n2c),(l2n2c),(c1dn2c),(c2n2c),(l1p2c),(l2p2c),(c1p2c),(c2p2c)]
# SNR [(lln2c),(rln2c),(cwdn2c),(llp2c),(rlp2c),(cwdp2c)]
header = ['hr', 'avail-p', 'inorg_n', 'inorg_p', 'sorbed_p', 'cl1',
          'cl2', 'cs1', 'cs2', 'snr1-cl1_n2c', 'snr2-cl2_n2c', 'snr3', 'snr4', 'snr5', 'snr6', 'snr7', 'snr8']


def est_uptake(snr, cl, cs, nt='n'):

    upt = 0
    if nt == 'n':
        data = snr[0:4]
    else:
        data = snr[4:]
    for x in range(4):
        if x < 2:
            upt += data[x] * cl[x]
        elif x > 1:
            upt += data[x] * cs[x - 2]
    return np.random.uniform(0, 2) * 1e-6 * upt


with open('carbon3_test.csv', 'w') as fh:
    CSV_WRITER = csv.writer(fh, delimiter=',')
    CSV_WRITER.writerow(header)
    # start inputs
    lnr = zeroes(6) + 0.005 * np.random.randint(0, 2)
    for x in range(8000):

        # nupt = np.random.uniform(0, 2) * 7.4e-3
        # pupt = np.random.uniform(0, 2) * 7.4e-3
        # UPDATE_ input VARIABLES
        snr_in = funcs.sp_snr
        avail_p = funcs.sp_available_p
        inorg_n = funcs.sp_available_n
        inorg_p = funcs.sp_in_p
        sorbed_p = funcs.sp_so_p
        cl = funcs.sp_csoil[0:2]
        cs = funcs.sp_csoil[2:]
        nupt = est_uptake(snr_in, cl, cs, 'n')
        pupt = est_uptake(snr_in, cl, cs, 'p')
        df = carbon3.carbon3(23, 0.8, 0, 0, 0, lnr, cl, cs,
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
        funcs.sp_available_n = df[1]  # - nupt
        funcs.sp_available_p = df[0]  # - pupt
        funcs.sp_in_p = df[2]
        funcs.sp_so_p = df[3]


data = pd.read_csv('carbon3_test.csv')
os.system('rm -rf carbon3_test.csv')
