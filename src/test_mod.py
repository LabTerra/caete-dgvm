
#   ____    _    _____ _____ _/\__
#  / ___|  / \  | ____|_   _| ____|
# | |     / _ \ |  _|   | | |  _|
# | |___ / ___ \| |___  | | | |___
#  \____/_/   \_\_____| |_| |_____|

import csv
import pickle
import glob
# import multiprocessing as mp
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from caete_input import data_in
# our module
import caete
from caete_module import global_par as gp
import os

# Example of (x, y) gridcell position
# ISIMIP2b input data
# manaus = (239, 183)
# hl = (533, 51)


# Creation of gridcell instances
grda = caete.gridcell_dyn(239, 183)
# grdb = caete.gridcell_dyn(219, 184)
# grdc = caete.gridcell_dyn(404, 186)
# grd1 = caete.gridcell_dyn(440, 61)
# grd2 = caete.gridcell_dyn(572, 55)

# Open input data
# The data_in-instance.pkl file is a instance of caete_inputs.py data_in saved as a pickle file
with open('../input/data_in-instance.pkl', 'rb') as fh:
    dt1 = pickle.load(fh)

# TODO implent dictionary creation: OK: Implemented in the create_input.py file

# Update(fill with input data) the created gridcells - init the gridcell with INPUT DATA
grda.init_caete_dyn(dt1, 'AM - Manaus')
# grdb.init_caete_dyn(dt1, 'AM - Oeste')
# grdc.init_caete_dyn(dt1, 'Africa - SAV-FOR')
# grd1.init_caete_dyn(dt1, 'RU - Oeste')
# grd2.init_caete_dyn(dt1, 'RU - Norte')

# A list of gridcells prepared to be simulated
# works = [grda, grdb, grdc, grd1, grd2]  # , grd3]


def f0(grd):
    """ Helper function. To be mapped over works"""
    # Apply the run_dyn function
    caete.run_dyn(grd)

    return grd


# Serial procesing: only one gridcell
grd = f0(grda)


# Parallel processing
# if __name__ == "__main__":

#     import multiprocessing as mp

# with mp.Pool(processes=5) as p:
#     result = p.map(f0, works)

fname = glob.glob1(os.getcwd(), "pools_din*")[0]
l = []
with open(fname, mode='r') as fh:
    for row in fh.readlines():
        l.append([x for x in row.rstrip().lstrip().split(' ') if x != ''])


def wcsv(li, mon=True):
    with open('pools.csv', 'w', newline='') as fh:
        wt = csv.writer(fh, delimiter=',')
        if mon:
            wt.writerow(['n', 'month', 'leaf', 'wood', 'root', 'litter', 'soil', 'wsoil',
                         'gpp', 'ra', 'npp', 'hr', 'rm', 'rg', 'wue', 'cue'])
        else:
            wt.writerow(['day', 'leaf', 'wood', 'root', 'litter1', 'litter2', 'soil1',
                         'soil2', 'wsoil', 'gpp', 'ra', 'npp', 'hr', 'rm', 'rg', 'wue', 'cue', 'rcm', 'ls'])

        wt.writerows(li)


s
wcsv(l, False)

data = pd.read_csv('pools.csv')


# --------------------------------------------------
# header = ['run', 'cleaf',
#           'cawood_comm',
#           'cfroot_comm',
#           'c_litter1',
#           'c_litter2',
#           'c_soil1',
#           'c_soil2',
#           'wsoil_comm',
#           'photo_comm',
#           'aresp_comm',
#           'npp_comm',
#           'het_resp',
#           'rm_comm',
#           'rg_comm',
#           'wue',
#           'cue',
#           'rcm_comm',
#           'life_strategy']

# # # def plots(): simple plot with matplotlib
# colors = ['g', 'r', 'b', 'm', 'y', 'k']
# legend = []

# for y in range(5):
#     X = []
#     l = []
#     legend.append(result[y].name)
#     for x in range(5):
#         if x == 0:
#             l.append(caete.gp.npls)
#         else:
#             l.append((result[y].area[x] > 0).sum())

#         X.append(x * caete.gp.nt1)
#     plt.plot(X, l, '-%s' % colors[y])
# plt.legend(legend)
# plt.xlabel('Day')
# plt.ylabel(' Número de Estratégias de Vida')
# plt.savefig("pls_surviving.png")

# k + run &
#                  & ,cleaf_cwm(k)&
#                  & ,cawood_cwm(k)&
#                  & ,cfroot_cwm(k)&
#                  & ,csoil_in(1)&
#                  & ,csoil_in(2)&
#                  & ,csoil_in(3)&
#                  & ,csoil_in(4)&
#                  & ,wsoil_cwm(k)&
#                  & ,photo_cwm(k)&
#                  & ,aresp_cwm(k)&
#                  & ,npp_cwm(k)&
#                  & ,hr_cwm(k)&
#                  & ,rm_cwm(k)&
#                  & ,rg_cwm(k)&
#                  & ,wue(k)&
#                  & ,cue(k)&
#                  & ,rcm_cwm(k)&
#                  & ,ls
