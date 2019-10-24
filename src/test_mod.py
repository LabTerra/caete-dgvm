
# import csv
import pickle
# ls 
# import multiprocessing as mp
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from caete_input import data_in
# our module
import caete
from caete_module import global_par as gp

# Example of (x, y) gridcell position
# ISIMIP2b input data
# manaus = (239, 183)
# hl = (533, 51)



# Creation of gridcell instances 
grda = caete.gridcell_dyn(239, 183)
grdb = caete.gridcell_dyn(219, 184)
grdc = caete.gridcell_dyn(404, 186) # aFRICA
grd1 = caete.gridcell_dyn(440, 61)
grd2 = caete.gridcell_dyn(572, 55)
# grd3 = caete.gridcell_dyn(239, 183)

# Open input data
# The test_data.pkl file is a instance of caete.data_in saved as a pickle file
with open('../input/data_in-instance.pkl', 'rb') as fh:
    dt1 = pickle.load(fh)

# TODO implent dictionary creation: OK: Implemented in the create_input.py file

# Update(fill with input data) the created gridcells - init the gridcell with INPUT DATA
grda.init_caete_dyn(dt1, 'AM - Manaus')
grdb.init_caete_dyn(dt1, 'AM - Oeste')
grdc.init_caete_dyn(dt1, 'Africa - SAV-FOR')
grd1.init_caete_dyn(dt1, 'RU - Oeste') 
grd2.init_caete_dyn(dt1, 'RU - Norte')


# caete.init_caete_dyn(grd3, dt1)
# caete.init_caete_dyn(grda, dt1, 'AM - Manaus')
# caete.init_caete_dyn(grdb, dt1, 'AM - Oeste')
# caete.init_caete_dyn(grdc, dt1, 'Africa - SAV-FOR')
# caete.init_caete_dyn(grd1, dt1, 'RU - Oeste') 
# caete.init_caete_dyn(grd2, dt1, 'RU - Norte')



## Make some expriments. Maybe
# grd1.tas = grdb.tas + (grdb.tas * 0.2)
# grd1.pr = grdb.pr - (grdb.pr * 0.80)
# grd2.tas = grdb.tas  +5.0
# grd2.pr = grdb.pr + (grdb.pr * 0.80) 
# grd3.tas = grdb.tas - 12
# grd3.pr = grdb.pr + 5e-5
# grdb.tas = grdb.tas - 2.0
# grdc.pr = grdc.pr - (grdb.pr * 0.50)
# grdc.tas = grdc.tas + 5
# grdc.pr = grdc.pr + 1e-4

# A list of gridcells prepared to be simulated
works = [grda, grdb, grdc, grd1, grd2] #, grd3]

def f0(grd):
    
    """ Helper function. To be mapped over works"""
    # Apply the run_dyn function
    caete.run_dyn(grd)

    return grd


# Example of gridcell execution:
# grd = f0(grda)

### Parallel processing 
#TODO
# Implement here gridcell communication?

# with mp.Pool(processes=5) as p:
#     result = p.map(f0, works)


# # A simple plot with matplotlib.
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
#     plt.plot(X, l, '-%s'%colors[y])
# plt.legend(legend)
# plt.xlabel('Day')
# plt.ylabel(' Número de Estratégias de Vida') 
# plt.savefig("pls_surviving.png")
