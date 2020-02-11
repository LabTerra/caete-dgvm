import sys
sys.path.insert(0, '../src/')
import numpy as np
import matplotlib.pyplot as plt
import caete_module as m
import plsgen as pls

# real(kind=r_4),dimension(ntraits,npft),intent(in) :: dt
dt = pls.table_gen(m.global_par.npls)

n_min = 0.3775999
p_min = 0.243

sto = np.zeros(shape=(3,m.global_par.npls), order='F')
# !OLD
# real(kind=r_8),dimension(npft),intent(inout) :: cl1_pft  ! initial BIOMASS cleaf compartment
cl1 = np.zeros(m.global_par.npls) + 1 * np.random.random(m.global_par.npls,) # np.random.random(m.global_par.npls) + 0.7
cf1 = np.zeros(m.global_par.npls) + 0.5 * np.random.random(m.global_par.npls,)#np.random.random(m.global_par.npls) + 0.7
ca1 = np.zeros(m.global_par.npls) + 0.1 * np.random.random(m.global_par.npls,)
# real(kind=r_8),dimension(npft),intent(inout) :: cf1_pft  !                 froot
# real(kind=r_8),dimension(npft),intent(inout) :: ca1_pft  !                 cawood
dleaf = np.random.random(m.global_par.npls) * 0.002
droot = np.random.random(m.global_par.npls) * 0.002
dwood = np.random.random(m.global_par.npls) * 0.002

io0  = [sto,cl1,ca1,cf1,dleaf,droot,dwood]
# nmin = np.array([0.5,])
#     # real(kind=r_8),intent(in) :: labile_p
# plab = np.array([0.5,])

def testb():
    # global nmin
    # global plab
    # real(kind=r_4),dimension(npft),intent(in) :: w1   !Initial (previous month last day) soil moisture storage (mm)
    w1 = np.zeros(m.global_par.npls) + 20
    # real(kind=r_4),dimension(npft),intent(in) :: g1   !Initial soil ice storage (mm)
    ga1 = np.zeros(m.global_par.npls)
    # real(kind=r_4),dimension(npft),intent(in) :: s1   !Initial overland snow storage (mm)
    s1 = np.zeros(m.global_par.npls)
    # real(kind=r_4),intent(in) :: ts                   ! Soil temperature (oC)
    tsoil = 26
    # real(kind=r_4),intent(in) :: temp                 ! Surface air temperature (oC)
    temp = 22
    # real(kind=r_4),intent(in) :: prec                 ! Precipitation (mm/monyh)
    prec = 5
    # real(kind=r_4),intent(in) :: p0                   ! Surface pressure (mb)
    p0 = 1013.25
    # real(kind=r_4),intent(in) :: ipar                 ! Incident photosynthetic active radiation mol Photons m-2 s-1
    ipar = 290 / 2.18e5
    # real(kind=r_4),intent(in) :: rh                   ! Relative humidity
    rh = 0.80
    # ! State variables INPUTS & OUTPUTS
    # !NEW
    # real(kind=r_8),intent(in) :: mineral_n   ! Mineral pools (SOIL) kg(Nutrient) m-2

    # real(kind=r_8),dimension(npft,3),intent(inout)  :: sto_budg ! Rapid Storage Pool (C,N,P)
    # real(kind=r_8),dimension(npft),intent(inout) :: dleaf  ! CHANGE IN cVEG (DAILY BASIS) TO GROWTH RESP
    # real(kind=r_8),dimension(npft),intent(inout) :: droot
    # real(kind=r_8),dimension(npft),intent(inout) :: dwood

    # subroutine daily_budget(dt,w1,g1,s1,ts,temp,prec,p0,ipar,rh&
    #                       &,mineral_n,labile_p,sto_budg&
    #                       &,cl1_pft,ca1_pft,cf1_pft,dleaf,dwood

    out = m.budget.daily_budget(dt,w1,ga1,s1,tsoil,temp,prec,p0,ipar,rh,sto,cl1,ca1,cf1,dleaf,droot,dwood)
    return out


def test_b_dyn():
    sto = np.zeros(shape=(3,m.global_par.npls), order='F')
    # !OLD
    # real(kind=r_8),dimension(npft),intent(inout) :: cl1_pft  ! initial BIOMASS cleaf compartment
    cl1 = np.zeros(m.global_par.npls) + 1 # np.random.random(m.global_par.npls) + 0.7
    cf1 = np.zeros(m.global_par.npls) + 1 #np.random.random(m.global_par.npls) + 0.7
    ca1 = np.zeros(m.global_par.npls) + 1
    # real(kind=r_8),dimension(npft),intent(inout) :: cf1_pft  !                 froot
    # real(kind=r_8),dimension(npft),intent(inout) :: ca1_pft  !                 cawood
    dleaf = np.random.random(m.global_par.npls) * 0.002
    droot = np.random.random(m.global_par.npls) * 0.002
    dwood = np.random.random(m.global_par.npls) * 0.002

    io0  = [sto,cl1,ca1,cf1,dleaf,droot,dwood]

    # out = m.budget.daily_budget(dt,w1,ga1,s1,tsoil,temp,prec,p0,ipar,rh,nmin,plab,sto,cl1,ca1,cf1,dleaf,droot,dwood)
    data = dict()
    w1 = np.zeros(m.global_par.npls) + 0.01
    # real(kind=r_4),dimension(npft),intent(in) :: g1   !Initial soil ice storage (mm)
    ga1 = np.zeros(m.global_par.npls)
    # real(kind=r_4),dimension(npft),intent(in) :: s1   !Initial overland snow storage (mm)
    s1 = np.zeros(m.global_par.npls)
    # real(kind=r_4),intent(in) :: ts                   ! Soil temperature (oC)
    tsoil = 22
    # real(kind=r_4),intent(in) :: temp                 ! Surface air temperature (oC)
    temp = 22
    # real(kind=r_4),intent(in) :: prec                 ! Precipitation (mm/monyh)
    prec = 203
    # real(kind=r_4),intent(in) :: p0                   ! Surface pressure (mb)
    p0 = 1013.25
    # real(kind=r_4),intent(in) :: ipar                 ! Incident photosynthetic active radiation mol Photons m-2 s-1
    ipar = 290 / 2.18e5
    # real(kind=r_4),intent(in) :: rh                   ! Relative humidity
    rh = 0.70

    out = m.budget.daily_budget(dt,w1,ga1,s1,tsoil,temp,prec,p0,ipar,rh,sto,cl1,ca1,cf1,dleaf,droot,dwood)
    ion  = [sto,cl1,ca1,cf1,dleaf,droot,dwood]
    io0  = [sto[:],cl1[:],ca1[:],cf1[:],dleaf[:],droot[:],dwood[:]]



    for nz in range(50):


        out = m.budget.daily_budget(dt,w1,ga1,s1,tsoil,temp,prec,p0,ipar,
            rh,io0[0],io0[1],io0[2],io0[3],io0[4],io0[5],io0[6],io0[7],io0[8])

        io0 = [io0[0],io0[1],io0[2],io0[3],io0[4],io0[5],io0[6],io0[7],io0[8]]

        #print(io0)
        inouts = 'sto,cl1,ca1,cf1,dleaf,droot,dwood'.split(',')

        lst = ('w2,g2,s2,smavg,ruavg,evavg,epavg,phavg,aravg,nppavg,laiavg,rcavg,f5avg,rmavg,\
               rgavg,cleafavg_pft,cawoodavg_pft,cfrootavg_pft,ocpavg,wueavg,cueavg,c_defavg,vcmax,\
               specific_la,nupt,pupt,litter_l,cwd,litter_fr,lnr')

        lst = lst.split(',')

        for x in range(9):
            pass
            print(inouts[x],' -> ',io0[x])
        for i,x in enumerate(out):
            pass
            print(lst[i], ' -> ', x)

        data[str(nz)] = (io0,out)

    return data



data = testb()
io = [sto,cl1,ca1,cf1,dleaf,droot,dwood]
#print(io0)
inouts = 'sto,cl1,ca1,cf1,dleaf,droot,dwood'.split(',')


lst = ('w2,g2,s2,smavg,ruavg,evavg,epavg,phavg,aravg,nppavg,laiavg,rcavg,f5avg,rmavg,rgavg,cleafavg_pft,cawoodavg_pft,cfrootavg_pft,ocpavg,wueavg,cueavg,c_defavg,vcmax,specific_la,nupt,pupt,litter_l,cwd,litter_fr,lnr')
lst = lst.split(',')

for x in range(7):
    print(inouts[x],' -> ',io[x])

for i,x in enumerate(data):
    print(lst[i], ' -> ', x)


data = dict(zip(lst,data))

