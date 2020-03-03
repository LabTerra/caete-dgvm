"""
Copyright 2017-2018 LabTerra

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""

# This module contains the procedures to create the set of PLant life strategies for CAETÊ runs

import os
from math import ceil
import csv
import numpy as np
from caete_module import photo as model

# Cache data of allocation coefficients
woody_allocations_file = "wallo.npy"
grassy_allocations_file = "gallo.npy"


def vec_ranging(values, new_min, new_max):
    """ Ranges the vector (1D) values(np.array) to min max
        - Normalize values
    """

    output = []
    old_min, old_max = min(values), max(values)

    for v in values:
        new_v = (new_max - new_min) / (old_max -
                                       old_min) * (v - old_min) + new_min
        output.append(new_v)

    return np.array(output, dtype=np.float32)


def check_viability(trait_values, wood):
    """ Check the viability of allocation(a) & residence time(ŧ) combinations.
        Some PLS combinations of allocation coefficients and residence times
        are not 'biomass acumulators' at low npp (< 0.01 kg m⁻² year⁻¹)
        do not have enough mass of carbon (< 0.01 kg m⁻²) in all CVEG compartments

        trait_values: np.array(shape=(6,), dtype=f32) allocation and residence time combination (possible PLS)
        wood: bool  Is this a woody PLS?
    """

    assert wood is not None
    rtur = np.array(model.spinup3(0.01, trait_values))
    if wood:
        if rtur[0] <= 0.01 or rtur[1] <= 0.01 or rtur[2] <= 0.01:
            return False
        return True
    else:
        if rtur[0] <= 0.01 or rtur[1] <= 0.01:
            return False
        return True


def assertion_data_size(dsize):
    """ Assertion of datasets sizes """

    g2w_ratio = 0.07
    diffg = ceil(dsize * g2w_ratio)
    diffw = int(dsize - diffg)
    assert diffg + diffw == dsize
    return diffg, diffw


def turnover_combinations(verbose=False):
    """CREATE the residence time and allocation combinations"""

    # constrained distributions (must sum up to 1.)

    if os.path.exists(grassy_allocations_file):

        plsa_grass = np.load(grassy_allocations_file)

    else:

        aleafg = np.arange(15., 86.0, 0.2000)
        arootg = np.arange(15., 86.0, 0.2000)

        plsa_grass = [[a / 100.0, 0.0, c / 100.0] for a in aleafg
                      for c in arootg if abs((a + c) - 100.0) < 0.00001]

        np.save(grassy_allocations_file, plsa_grass)

    if os.path.exists(woody_allocations_file):

        plsa_wood = np.load(woody_allocations_file)

    else:

        aleafw = np.arange(15., 86.0, 0.2000)
        arootw = np.arange(15., 86.0, 0.2000)
        awood = np.arange(15., 86.0, 0.2000)

        plsa_wood = [[a / 100.0, b / 100.0, c / 100.0] for a in aleafw for b in awood
                     for c in arootw if abs((a + b + c) - 100.) < 0.00001]
        np.save(woody_allocations_file, plsa_wood)

    if verbose:
        print('Number of combinations = %d' %
              (len(plsa_grass) + len(plsa_wood)))

    return np.array(plsa_wood), np.array(plsa_grass)


def table_gen(NPLS):
    """AKA main - generate a trait table for CAETÊ - save it to a .csv"""

    diffg, diffw = assertion_data_size(NPLS)
    plsa_wood, plsa_grass = turnover_combinations(True)

    alloc_w = []
    alloc_g = []
    r_ceil = 30000

# REVER O TEMPO DE RESIDÊNCIA DAS RAÌZES FINAS - VARIAR ENTRE 1 mes e 2 anos
    index0 = 0
    rtime = vec_ranging(np.random.normal(0.0, 1.0, r_ceil), 0.08333333334, 8.3)
    while index0 < diffg:
        restime = np.zeros(shape=(3,), dtype=np.float32)

        allocatio = plsa_grass[np.random.randint(0, plsa_grass.shape[0])]
        restime[0] = rtime[np.random.randint(0, r_ceil)]
        restime[1] = 0.0
        restime[2] = rtime[np.random.randint(0, r_ceil)]

        data_to_test0 = np.concatenate((restime, allocatio), axis=0,)
        if check_viability(data_to_test0, False):
            alloc_g.append(data_to_test0)
            index0 += 1
    # Creating woody plants (maybe herbaceous)
    index1 = 0
    rtime_wood = vec_ranging(np.random.normal(0, 1, r_ceil), 1., 100.)
    while index1 < diffw:
        restime = np.zeros(shape=(3,), dtype=np.float32)
        allocatio = plsa_wood[np.random.randint(0, plsa_wood.shape[0])]
        restime[0] = rtime[np.random.randint(0, r_ceil)]
        restime[1] = rtime_wood[np.random.randint(0, r_ceil)]
        restime[2] = rtime[np.random.randint(0, r_ceil)]
        data_to_test1 = np.concatenate((restime, allocatio), axis=0,)
        if check_viability(data_to_test1, True):
            alloc_w.append(data_to_test1)
            index1 += 1

    alloc_g = np.array(alloc_g)
    alloc_w = np.array(alloc_w)

    alloc = np.concatenate((alloc_g, alloc_w), axis=0,)

    # # # COMBINATIONS
    # # # Random samples from  distributions (g1, tleaf ...)
    # # # Random variables
    g1 = np.random.uniform(1.0, 15.0, NPLS)
    # g1 = vec_ranging(np.random.beta(1.2, 2, NPLS), 1.0, 15.0) # dimensionles
    # # vcmax = np.random.uniform(3e-5, 100e-5,N) # molCO2 m-2 s-1
    vcmax = np.zeros(NPLS,)

    # # C4 STYLE
    c4 = np.zeros((NPLS,), dtype=np.float32)
    n123 = ceil(alloc_g.shape[0] * 0.70)
    c4[0:n123 - 1] = 1.0

    # # Nitrogen and Phosphorus content in carbon pools
    # # C : N : P
    # # ----
    leaf_n2c = np.random.uniform(0.005, 0.01, NPLS)
    awood_n2c = np.random.uniform(0.005, 0.01, NPLS)
    awood_n2c[0:alloc_g.shape[0]] = 0.0
    froot_n2c = np.random.uniform(0.005, 0.01, NPLS)

    leaf_p2c = np.random.uniform(0.001, 0.01, NPLS)
    awood_p2c = np.random.uniform(0.001, 0.01, NPLS)
    awood_p2c[0:alloc_g.shape[0]] = 0.0
    froot_p2c = np.random.uniform(0.001, 0.01, NPLS)

    stack = (g1, vcmax, alloc[:, 0], alloc[:, 1], alloc[:, 2],
             alloc[:, 3], alloc[:, 4], alloc[:, 5], c4, leaf_n2c,
             awood_n2c, froot_n2c, leaf_p2c, awood_p2c, froot_p2c)

    head = ['g1', 'vcmax', 'tleaf', 'twood', 'troot', 'aleaf', 'awood', 'aroot', 'c4',
            'leaf_n2c', 'awood_n2c', 'froot_n2c', 'leaf_p2c', 'awood_p2c', 'froot_p2c']
    # # for i,x in enumerate(head):
    # #     print(i + 1, ': ',x)

    pls_table = np.vstack(stack)
    # # pls_table_F = np.empty(shape=(15,n),order='F')
    # # pls_table_F = pls_table

    # # ___side_effects
    with open('pls_attrs.csv', mode='w') as fh:
        writer = csv.writer(fh, delimiter=',')
        writer.writerow(head)
        for x in range(pls_table.shape[1]):
            writer.writerow(list(pls_table[:, x]))
        # writer.writerows(pls_table)

    out_arr = np.asfortranarray(pls_table, dtype=np.float32)
    np.savetxt('pls_ex.txt', out_arr.T, fmt='%.12f')

    return out_arr
