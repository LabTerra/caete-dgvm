import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.insert(0, '../src/')

import plsgen as pls
import caete_module as m

dt = pls.table_gen(m.global_par.npls)


def print_alloc(trait):

    # real(kind=r_8),dimension(3),intent(out) :: storage_out
    # real(kind=r_8),intent(out) :: scl2 ! final carbon content on leaf compartment (KgC/m2)
    # real(kind=r_8),intent(out) :: sca2 ! final carbon content on aboveground woody biomass compartment (KgC/m2)
    # real(kind=r_8),intent(out) :: scf2 ! final carbon content on fine roots compartment (KgC/m2)
    # real(kind=r_8),intent(out) :: cwd  ! coarse wood debris (to litter)(C) g m-2
    # real(kind=r_8),intent(out) :: root_litter ! to litter g(C) m-2
    # real(kind=r_8),intent(out) :: leaf_litter ! to litter g(C) m-2
    # real(kind=r_8),intent(out) :: nout ! N plant uptake g(N) m-2
    # real(kind=r_8),intent(out) :: pout ! P plant uptake g(P) m-2
    # real(kind=r_8),dimension(6),intent(out) :: litter_nutrient_ratios

    nmin = 5.65e-3  # Kg m-2
    plab = 1.45e-3  # kg m-2
    l_init = 1  # Kg m-2
    r_init = 1  # Kg m-2
    w_init = 6  # Kg m-2

    if not dt[6][trait] > 0.0:
        w_init = 0.0

    sto = np.array([0.0, 0.0, 0.0])  # g m-2

    d = m.photo.allocation(dt[:, trait], 2.0, nmin,
                           plab, l_init, w_init, r_init, sto)

    variables = ['storage_pool', 'cleaf', 'cawood', 'cfroot', 'leaf_litter',
                 'cwd', 'root_litter', 'nup', 'pup', 'litter_nutrient_ratios']

    units = ['(g/m2)', 'KgC/m2', 'KgC/m2', 'KgC/m2', 'g(C)/m2', 'g(C)/m2',
             'g(C)/m2', 'g(N)/m2', 'g(P)/m2', 'g(Nutrient)/g(C)']

    for i, v in enumerate(variables):
        print(v, '--->', d[i], '-unit-->', units[i])
    print(' ![(lln2c),(rln2c),(cwdn2c),(llp2c),(rlp2c),(cwdp2c)]')


def test12():
    for x in range(1, m.global_par.npls):
        print_alloc(x)
