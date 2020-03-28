
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


experiment = ['inc_uptake', 'inc_litter_flow', 'dec_uptake', ]

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
header = ['Het. Respiration', 'Available N', 'Inorganic P',
          'Sorbed P', 'Available P', 'Litter I)', 'Litter II)', 'Soil I)', 'Soil II)',
          'N:C Litter I', 'N:C Litter II', 'N:C Soil I', 'N:C Soil II', 'P:C Litter I',
          'P:C Litter II', 'P:C Soil I', 'P:C Soil II']


def reset_soil():

    funcs.sp_snr = zeroes(8)
    funcs.sp_csoil[0:2] = zeroes(2) + 500
    funcs.sp_csoil[2:] = zeroes(2) + 500
    funcs.sp_available_n = 0.50
    funcs.sp_available_p = 0.20
    funcs.sp_in_p = 0.0
    funcs.sp_so_p = 0.0


reset_soil()


def est_uptake(lnr, ll, rl, cwd, nt='n'):

    upt = 0
    if nt == 'n':
        data = lnr[0:3]
        a = 1
    else:
        data = lnr[3:]
        a = 0.8
    for x in range(3):
        upt += data[x] * ll * a
        upt += data[x] * rl * a
        upt += data[x] * cwd * a
    return np.random.uniform(0, 2) * upt


with open('carbon3_test.csv', 'w') as fh:
    CSV_WRITER = csv.writer(fh, delimiter=',')
    CSV_WRITER.writerow(header)
    # start inputs
    # lnr = zeroes(6) + 0.01 * np.random.randint(0, 2)

    for x in range(40000):
        lnr_std = zeroes(6) + 0.03  # * np.random.randint(0, 2)
        lnr_std[3:] = lnr_std[3:] * 0.8
        if x > 50000:
            leaf_litter = np.random.uniform(0, 2) * 10 * 2
            root_litter = np.random.uniform(0, 2) * 10 * 2
            cwd = np.random.uniform(0, 2) * 10 * 2
        else:
            leaf_litter = np.random.uniform(0, 2) * 10
            root_litter = np.random.uniform(0, 2) * 10
            cwd = np.random.uniform(0, 2) * 10
        if x > 20000:
            lnr = zeroes(6) + 0.03  # * np.random.randint(0, 2)
            lnr[3:] = lnr[3:] * 0.8
        else:
            lnr = zeroes(6) + 0.003  # * np.random.randint(0, 2)
            lnr[3:] = lnr[3:] * 0.8

        snr_in = funcs.sp_snr
        avail_p = funcs.sp_available_p
        inorg_n = funcs.sp_available_n
        inorg_p = funcs.sp_in_p
        sorbed_p = funcs.sp_so_p
        cl = funcs.sp_csoil[0:2]
        cs = funcs.sp_csoil[2:]

        nupt = est_uptake(lnr_std, leaf_litter, root_litter, cwd, 'n')
        pupt = est_uptake(lnr_std, leaf_litter, root_litter, cwd, 'p')

        if x > 50000:
            pupt = pupt * 2
            nupt = nupt * 2

        df = carbon3.carbon3(23, 0.8, leaf_litter, root_litter, cwd, lnr, cl, cs,
                             snr_in, avail_p, inorg_n, inorg_p, sorbed_p)
        # print(df)
        line = [df[-1],
                df[1],
                df[2],
                df[3],
                df[0],
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
        funcs.sp_available_n = df[1] - nupt
        funcs.sp_available_p = df[0] - pupt
        funcs.sp_in_p = df[2]
        funcs.sp_so_p = df[3]


data = pd.read_csv('carbon3_test.csv')
os.system('rm -rf carbon3_test.csv')

dt = data.iloc(1)

n1 = 'n2c_inc_litter_Q.png'
n2 = 'p2c_inc_litter_Q.png'
n3 = 'in_nut_inc_litter_Q.png'
n4 = 'cpools_inc_litter_Q.png'
n5 = 'in_nut_comp_inc_litter_Q.png'

dt[9:13].plot(cmap=plt.get_cmap('cool'))
plt.xlabel('Iterations')
plt.ylabel('g(N) g(C)⁻¹')
plt.savefig(n1, dpi=700)
plt.show()

dt[13:].plot(cmap=plt.get_cmap('rainbow'))
plt.xlabel('Iterations')
plt.ylabel('g(P) g(C)⁻¹')
plt.savefig(n2, dpi=700)
plt.show()

dt[1:5].rolling(1000).mean().plot(cmap=plt.get_cmap('cividis'))
plt.xlabel('Iterations')
plt.ylabel('g(Nutrient) m⁻²')
plt.savefig(n3, dpi=700)
plt.show()

dt[5:9:].plot(cmap=plt.get_cmap('cividis'))
plt.xlabel('Iterations')
plt.ylabel('g(C) m⁻²')
plt.savefig(n4, dpi=700)
plt.show()

dt[1:5].rolling(1000).mean().plot(
    subplots=True, figsize=(6, 6), cmap=plt.get_cmap('viridis'))
plt.savefig(n5, dpi=700)
plt.show()
