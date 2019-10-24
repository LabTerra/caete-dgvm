import numpy as np
import matplotlib.pyplot as plt
import caete_module as model
from math import exp
from math import log
from math import sqrt
t = True
f = False
c3 = np.zeros(2000)
c3_1 = np.zeros(2000)
c4 = np.zeros(2000)
c4_1 = np.zeros(2000)
c3a = np.zeros(2000)
c3a_1 = np.zeros(2000)
c4a = np.zeros(2000)
c4a_1 = np.zeros(2000)

temp = np.linspace(-10,50,2000)
f1 = model.water.available_energy

nitrogen = 0.1 # g m-2
phosphorus =  0.02 # g m-2

A = model.photo.photosynthesis_rate
#  model.water.available_energy(temp[x])
def plotcomp():
    nbio = nitrogen  
    pbio = phosphorus
    sto = np.array([0.0002, 0.000003, 0.00000018])
    for x in list(range(temp.size)):
        c3a[x], c3a_1[x], sdf1 = A(temp[x],993.25,f1(temp),f,0,nbio,pbio,1,sto)
        c4a[x], c4a_1[x], sdf2 = A(temp[x],993.25,f1(temp),f,1,nbio,pbio,1,sto)
    plt.plot(temp, c3a, 'g', temp, c4a, 'r')
    plt.legend(['C3','C4'])
    plt.xlabel(r'T (°C)')
    plt.ylabel(r'A (mol/m²/s)')
    plt.show()
    return sdf1, sdf2

  # subroutine photosynthesis_rate(temp,p0,ipar,ll,c3,nbio,pbio,&
                                # & leaf_turnover,f1ab,vm)

def a():
    pass

def test_nu_vcmax(nbio2, pbio2):
    # vm_nutri = 3.946 + 0.921 * log(nbio) - 0.121 * log(pbio)
    # vm_nutri = vm_nutri + 0.282 * log(nbio) * log(pbio)
    # vm_nutri = exp(vm_nutri)
    # # return vm_nutri 
    # vm_nutri = 3.946 + 0.921 * nbio - 0.121 * pbio
    # vm_nutri = vm_nutri + 0.282 * nbio * pbio
    # return vm_nutri * 1e-6

    # vm_nutri = 3.946 + (0.921 * nbio2) - (0.121 * pbio2)
    # vm_nutri = (vm_nutri + (0.282 * nbio2) * pbio2) * 1e-6
    # return vm_nutri #Vcmax

    vm_nutri = 3.946 + 0.921 * log(nbio2) - 0.121 * log(pbio2)
    vm_nutri = vm_nutri + 0.282 * log(nbio2) * log(pbio2)
    
    # vm = exp(vm_nutri) * 1e-6 ! Vcmax

    # vm_nutri = 3.946 + (0.921 * log(nbio2)) - (0.121 * log(pbio2))
    # vm_nutri = exp((vm_nutri + (0.282 * log(nbio2)) * log(pbio2)))
    
    # vm_nutri = 3.946 + (0.921 * log(nbio2)) - (0.121 * log(pbio2))
    # vm_nutri = (vm_nutri + (0.282 * log(nbio2)) * log(pbio2)) 
    vm_nutri = exp(vm_nutri) * 1e-6
    return vm_nutri


def test_c4(c1='g', c2='y'):
    sto = np.array([0.002, 0.0003, 0.00018])
    temp2 = 22.0
    p0 = 1013.25
    nbio = nitrogen
    pbio = phosphorus
    rad = np.linspace(0,700,2000) / 2.18e5
    x = 0
    for i in rad:
        c3[x], c3_1[x], sdf = model.photo.photosynthesis_rate(temp2,960,i,t,1,nbio,pbio,1.7,sto)
        c4[x], c4_1[x], sdf2 = model.photo.photosynthesis_rate(temp2,960,i,t,0,nbio,pbio,1.7,sto)
        x += 1
    plt.plot(rad, c3, c1, rad, c4, c2)
    plt.legend(['C3','C4'])
    plt.xlabel(r'PAR (mol/m²/s)')
    plt.ylabel(r'A (mol/m²/s)')
    plt.show()

def clean():
    global c3, c3_1
    global c4, c4_1

    c3, c3_1 = 0.0, 0.0
    c4, c4_1 = 0.0, 0.0

def vcmax_lim():
    vm = np.linspace(3e-5, 100e-5, 200)
    out = 6 * (1e4**(-1 * ((vm - 3e-5) / (100e-5 - 3e-5))) / 10.0)
    plt.plot(vm,out)


def c4photo_test(alphap,vm):
    """From Chen et al. 1994"""

    tk = 25 + 273.15  #! K
    t25 = 273.15 + 25   #! tk at 25°C  
    ipar1 = (250 / 2.18e5) * 1e6  #! µmol m-2 s-1 - 1e6 converts mol to µmol
    vpm25 =  85.69      #! µmol m-2 s-1 PEPcarboxylase CO2 saturated rate of carboxilation at 25°C
    h_vpm = 185075.0    #! Arrhenius eq. constant
    s_vpm = 591         #! Arrhenius eq. constant
    r_vpm = 8.314       #! Arrhenius eq. constant
    e_vpm = 60592.0     #! Arrhenius eq. constant
    #alphap = 0.0913     #! parameter for v4m  
    kp = 82.0           #! µmol mol-1


    dummy1 = 1 + exp((s_vpm * t25 - h_vpm)/(r_vpm * t25))
    dummy2 = 1 + exp((s_vpm * tk - h_vpm)/(r_vpm * tk))
    dummy0 = dummy1/dummy2 
    vpm =  vpm25 * exp((-e_vpm/r_vpm) * (1/tk - 1/t25)) * dummy0
       
    #       ! actual PEPcarboxylase rate 
    v4m = (alphap * ipar1) / sqrt(1 + alphap**2 * ipar1**2 / vpm**2)

    # mesophyl 
    cm0 = 1.674 - 6.1294 * 10**(-2) * temp
    cm1 = 1.1688 * 10**(-3) * temp ** 2
    cm2 = 8.8741 * 10**(-6) * temp ** 3
    cm = 0.7 * ca * ((cm0 + cm1 - cm2) / 0.73547)
       
    # ! Michaelis-Menten eq. PEPcarboxylase kinetics
    # ! We are assuming here that assimilation is the same of PEP carboxylation
    # ! I want to penalize this assimilation using the maximum Rubisco carboxilation Rate (vm).
    # ! we know that the maximum Vcmax allowed in the model is 100e-5 mol m-2 s-1. Minimum is 3e-5
    # ! Normalizing 
    # aux666 = 8 * (1e4**(-1 * ((vm - 3e-5) / (100e-5 - 3e-5))) / 10.0) # ! range 0 to 0.10 (aprox)  
    f1ab = real((V4m * cm) / (kp + cm),r_4) * 1e-6 

    return f1ab


def linreg():
    from scipy.stats import linregress
    #from pylab import plot, title, show , legend
    #Sample data
    vcmax = np.linspace(3e-5,25e-5,200)
    alphap = np.linspace(3e-5,25e-5,200)


    d = linregress(alphap,vcmax)
    # print('Linear regression using stats.linregress')
    # print('parameters: a=%.2f b=%.2f \nregression: a=%.2f b=%.2f, std error= %.3f' % (a,b,a_s,b_s,stderr))
    # #Linear regressison -polyfit - polyfit can be used other orders polys
    # (ar,br)=polyfit(vcmax,alphap,1)
    # xr=polyval([ar,br],t)
    # #compute the mean square error
    # err=sqrt(sum((xr-xn)**2)/n)

    # print('Linear regression using polyfit')
    # print('parameters: a=%.2f b=%.2f \nregression: a=%.2f b=%.2f, ms error= %.3f' % (a,b,ar,br,err))
    return d


def collatz_model():

    # A = Wp - L

    pass

def cm(temp,ca):

    cm0 = 1.674 - 6.1294 * 10**(-2) * temp
    cm1 = 1.1688 * 10**(-3) * temp ** 2
    cm2 = 8.8741 * 10**(-6) * temp ** 3
    cm = 0.5 * ca * ((cm0 + cm1 - cm2) / 1.2)
    cm0 = cm / 9.90111
    return cm, cm0

def nlignin(leaf_turnover):
    tl0 = model.photo.spec_leaf_area(1.0/12.0)
    tlm = model.photo.spec_leaf_area(8.3)
    # ! tlm --- 0.8
    # ! tl0 --- 0.2
    tl = model.photo.spec_leaf_area(leaf_turnover)

    tl = 1.0 - (((tl - tl0) / (tlm - tl0)) * (0.8 - 0.2) + 0.2)
    return tl

def sla(tau_leaf):
    leaf_t_months = tau_leaf*12. #! turnover time in months
    leaf_t_coeff = leaf_t_months/100. #!1 - 100 months == ~ 1/12 to 8.3 years (TRY-kattge et al. 2011; Jedi-Pavlick 2012) 
    # !<<<<<<< HEAD
    leaf_turnover =  (365.0/12.0) * exp(2.6*leaf_t_coeff)
    # !=======
    # !>>>>>>> 2ca26587106701e0883c85ea235c431dcb9ee97b
    # leaf_turnover =  (365.0/12.0) * (10 ** (2.0*leaf_t_coeff))
    # leaf_turnover =  (365.0/12.0) * (5.3 ** (leaf_t_coeff))
    return (3e-2 * (365.0/leaf_turnover)**(-1.02))

def srun():
    d = np.zeros(100,)
    t = np.linspace(1/12,8.3,100)

    for i,x in enumerate(t):
        # d[i] = sla(x)
        d[i] = model.photo.spec_leaf_area(x)

    plt.plot(t,d)
    plt.show()

def lai():
    pass
