# Function to test productivity
import sys

sys.path.insert(0, '../src/')

import numpy as np
import matplotlib.pyplot as plt
import caete_module as m
import plsgen as pls

dt = pls.table_gen(m.global_par.npls)


def test_prod(n):

    # real(kind=r_4), intent(in) :: temp                 !Mean monthly temperature (oC)
    temp = 22.0
    # real(kind=r_4), intent(in) :: p0                   !Mean surface pressure (hPa)
    p0 = 1010.0
    # real(kind=r_4), intent(in) :: w
    w = 50
    # real(kind=r_4), intent(in) :: ipar                 !Incident photosynthetic active radiation (w/m2)
    ipar = 160 / 2.18e5
    # real(kind=r_4), intent(in) :: rh,emax              !Relative humidity/MAXIMUM EVAPOTRANSPIRATION
    rh = 0.60
    emax = 1.5  # kg m-2 day-1
    # real(kind=r_8), intent(in) :: cl1, cf1, ca1        !Carbon in plant tissues (kg/m2)
    cl1 = np.zeros(3,) + 0.17
    cf1 = 0.5
    ca1 = 6
    # real(kind=r_8), intent(in) :: beta_leaf            !npp allocation to carbon pools (kg/m2/day)
    beta_leaf = 1e-2
    beta_awood = 1e-2
    beta_froot = 1e-2
    # real(kind=r_8), intent(in) :: beta_awood
    # real(kind=r_8), intent(in) :: beta_froot
    # real(kind=r_8), dimension(3), intent(in) :: sto1
    sto1 = np.array([1e-8, 0.000004, 0.00000001])
    # logical(kind=l_1), intent(in) :: light_limit                !True for no ligth limitation
    ll = True

    # IN
    # subroutine prod(dt,pft,light_limit,temp,p0,w,ipar,rh,emax,cl1,&
    #               & ca1,cf1,beta_leaf,beta_awood,beta_froot,sto1
    # OUT
    # ,ph,ar,nppa,laia,f5,f1,vpd,rm,rg,rc,wue,c_defcit,vm_out,sto2)
    output = m.productivity.prod(dt.T[n], ll, temp, p0, w, ipar, rh,
                                 emax, cl1, ca1, cf1, beta_leaf, beta_awood, beta_froot, sto1)

    d = 'ph,ar,nppa,laia,f5,vpd,rm,rg,rc,wue,c_defcit,vm_out,sla,sto2'.split(
        ',')

    for i, x in enumerate(d):
        print(x, ': ', output[i])
