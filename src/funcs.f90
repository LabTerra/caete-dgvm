! Copyright 2017- LabTerra

!     This program is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.

!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR 2PURPOSE.  See the
!     GNU General Public License for more details.

!     You should have received a copy of the GNU General Public License
!     along with this program.  If not, see <http://www.gnu.org/licenses/>.

module photo

   ! Module defining functions related with CO2 assimilation and other processes in CAETE
   ! Some of these functions are based in CPTEC-PVM2, others are new features

   implicit none
   private

   ! functions(f) and subroutines(s) defined here
   public ::                    &
        gross_ph               ,& ! (f), gross photosynthesis (kgC m-2 y-1)
        leaf_area_index        ,& ! (f), leaf area index(m2 m-2)
        f_four                 ,& ! (f), auxiliar function (calculates f4sun or f4shade or sunlai)
        spec_leaf_area         ,& ! (f), specific leaf area (m2 g-1)
        water_stress_modifier  ,& ! (f), F5 - water stress modifier (dimensionless)
        photosynthesis_rate    ,& ! (s), leaf level CO2 assimilation rate (molCO2 m-2 s-1)
        canopy_resistence      ,& ! (f), Canopy resistence (from Medlyn et al. 2011a) (s/m) == m s-1
        stomatal_conductance   ,& ! (f), IN DEVELOPMENT - return stomatal conductance
        vapor_p_defcit         ,& ! (f), Vapor pressure defcit  (kPa)
        tetens                 ,& ! (f), Maximum vapor pressure (hPa)
        nrubisco               ,& ! (f), Fraction of N not in lignin (disponible to rubisco)
        nlignin                ,& ! (f), Fraction of N in lignin
        m_resp                 ,& ! (f), maintenance respiration (plants)
        spinup2                ,& ! (s), SPINUP function for CVEG pools
        spinup3                ,& ! (s), SPINUP function to check the viability of Allocation/residence time combinations
        g_resp                 ,& ! (f), growth Respiration (kg m-2 yr-1)
        pft_area_frac          ,& ! (s), area fraction by biomass
        allocation             ,& ! (s), daily npp allocation subroutine
        water_ue               ,& ! (f), water use efficiency
        rc_canopy_scaling         ! (f), Try it!

contains

   !=================================================================
   !=================================================================

   function gross_ph(f1,cleaf,sla) result(ph)
      ! Returns gross photosynthesis rate (kgC m-2 y-1) (GPP)
      use types, only: r_4, r_8
      !implicit none

      real(r_4),intent(in) :: f1    !molCO2 m-2 s-1
      real(r_8),intent(in) :: cleaf !kgC m-2
      real(r_8),intent(in) :: sla   !m2 gC-1
      real(r_4) :: ph

      real(r_4) :: f4sun
      real(r_4) :: f4shade

      f4sun = f_four(1,cleaf,sla)
      f4shade = f_four(2,cleaf,sla)

      ph = 0.012*31557600.0*f1*f4sun*f4shade

   end function gross_ph

   !=================================================================
   !=================================================================

   function leaf_area_index(cleaf, sla) result(lai)
      ! Returns Leaf Area Index m2 m-2

      use types, only: r_4, r_8
      !implicit none

      real(r_8),intent(in) :: cleaf !kgC m-2
      real(r_8),intent(in) :: sla   !m2 gC-1
      real(r_8) :: lai


      lai  = cleaf * 1000. * sla  ! Converts cleaf from (KgC m-2) to (gCm-2)


   end function leaf_area_index

   !=================================================================
   !=================================================================

   function spec_leaf_area(tau_leaf) result(sla)
      ! based on JeDi DGVM
      use types, only : r_4, r_8
      !implicit none

      real(r_4),intent(in) :: tau_leaf  !years
      real(r_4) :: sla   !m2 gC-1

      real(r_4) :: leaf_t_months
      real(r_4) :: leaf_t_coeff
      real(r_4) :: leaf_turnover

      leaf_t_months = tau_leaf*12. ! turnover time in months
      leaf_t_coeff = leaf_t_months/100. !1 - 100 months == ~ 1/12 to 8.3 years (TRY-kattge et al. 2011; Jedi-Pavlick 2012)

      ! leaf_turnover =  (365.0/12.0) * exp(2.6*leaf_t_coeff)
      leaf_turnover =  (365.0/12.0) * (10 ** (2.0*leaf_t_coeff))

      ! sla = (3e-2 * (365.0/leaf_turnover)**(-1.02))
      sla = (3e-2 * (365.0/leaf_turnover)**(-0.46))

   end function spec_leaf_area

   !=================================================================
   !=================================================================

   function f_four(fs,cleaf,sla) result(lai_ss)
      ! Function used to scale LAI from leaf to canopy level (2 layers)
      use types, only: i_4, r_4, r_8
      use photo_par, only: p26, p27
      !implicit none

      integer(i_4),intent(in) :: fs !function mode:
      ! 1  == f4sun   --->  to gross assimilation
      ! 2  == f4shade --->  too
      ! 90 == sun LAI
      ! 20 == shade LAI
      ! Any other number returns sunlai (not scaled to canopy)

      real(r_8),intent(in) :: cleaf ! carbon in leaf (kg m-2)
      real(r_8),intent(in) :: sla   ! specific leaf area (m2 gC-1)
      real(r_4) :: lai_ss           ! leaf area index (m2 m-2)

      real(r_8) :: lai
      real(r_8) :: sunlai
      real(r_8) :: shadelai

      lai = leaf_area_index(cleaf,sla)

      sunlai = (1.0D0-(dexp(-p26*lai)))/p26
      shadelai = lai - sunlai

      lai_ss = real(sunlai,r_4)

      if (fs .eq. 90) then
         return
      endif
      if (fs .eq. 20) then
         lai_ss = real(shadelai, r_4)
         return
      endif

      !Scaling-up to canopy level (dimensionless)
      !------------------------------------------
      !Sun/Shade approach to canopy scaling !Based in de Pury & Farquhar (1997)
      !------------------------------------------------------------------------
      if(fs .eq. 1) then
         ! f4sun
         lai_ss = real((1.0-(dexp(-p26*sunlai)))/p26,r_4) !sun decl 90 degrees
         return
      endif

      if(fs .eq. 2) then
         !f4shade
         lai_ss = real((1.0-(dexp(-p27*shadelai)))/p27,r_4) !sun decl ~20 degrees
         return
      endif
   end function f_four

   !=================================================================
   !=================================================================

   function rc_canopy_scaling(rc) result(canopy_rc)
      use types
      use photo_par, only: p26, p27

      real(r_4), intent(in) :: rc
      real(r_4), dimension(2) :: canopy_rc
      real(r_4) :: rc_sun, rc_sha

      rc_sun = real((1.0-(exp(-p26*rc/500.0)))/p26, r_4)
      rc_sha = rc/500.0 - rc_sun

      canopy_rc(1) = real((1.0-(exp(-p26*rc_sun)))/p26,r_4)    ! frc_sun
      canopy_rc(2) = real((1.0-(exp(-p27*rc_sha)))/p27,r_4)    ! frc_sha
      return
   end function rc_canopy_scaling

   !=================================================================
   !=================================================================

   function water_stress_modifier(w, cfroot, rc, ep) result(f5)
      use types, only: r_4, r_8
      use global_par, only: csru, wmax, alfm, gm, rcmin
      !implicit none

      real(r_4),intent(in) :: w      !soil water mm
      real(r_8),intent(in) :: cfroot !carbon in fine roots kg m-2
      real(r_4),intent(in) :: rc     !Canopy resistence 1/(micromol(CO2) m-2 s-1)
      real(r_4),intent(in) :: ep     !potential evapotranspiration
      real(r_4) :: f5


      real(r_8) :: pt
      real(r_8) :: gc
      real(r_8) :: wa
      real(r_8) :: d
      real(r_8) :: f5_64

      wa = w/wmax

      pt = csru*(cfroot*1000.) * wa  !(based in Pavlick et al. 2013; *1000. converts kgC/m2 to gC/m2)
      if(rc .gt. 0.0) then
         gc = (1.0/(rc * 1.15741e-08))  ! s/m
      else
         gc =  1.0/(rcmin * 1.15741e-08) ! BIANCA E HELENA - Mudei este esquema..
      endif

      !d =(ep * alfm) / (1. + gm/gc) !(based in Gerten et al. 2004)
      d = (ep * alfm) / (1. + (gm/gc))
      if(d .gt. 0.0) then
         f5_64 = pt/d
         f5_64 = exp(-0.1 * f5_64)
         f5_64 = 1.0 - f5_64
      else
         f5_64 = wa
      endif

      f5 = real(f5_64,4)
      if (f5 .lt. 0.0) f5 = 0.0
   end function water_stress_modifier

   ! =============================================================
   ! =============================================================

   function canopy_resistence(vpd_in,f1_in,g1) result(rc2_in)
      ! return stomatal resistence based on Medlyn et al. 2011a
      ! Coded by Helena Alves do Prado

      use types, only: r_4 ,r_8
      use global_par, only: ca


      !implicit none

      real(r_4),intent(in) :: f1_in    !Photosynthesis (molCO2/m2/s)
      real(r_4),intent(in) :: vpd_in   !hPa
      real(r_4),intent(in) :: g1       ! model m (slope) (sqrt(kPa))
      real(r_4) :: rc2_in              !Canopy resistence (sm-1)

      !     Internal
      !     --------
      real(r_8) :: gs       !Canopy conductance (molCO2 m-2 s-1)
      real(r_8) :: D1       !sqrt(kPA)
      real(r_4) :: vapour_p_d

      vapour_p_d = vpd_in
      ! Assertions
      if(vpd_in .le. 0.0) vapour_p_d = 0.001
      if(vpd_in .gt. 4.0) vapour_p_d = 4.0
      ! print *, 'vpd going mad in canopy_resistence'
      ! stop
      ! endif

      D1 = sqrt(vapour_p_d)
      gs = 0.01 + 1.6 * (1.0 + (g1/D1)) * ((f1_in * 1e6)/ca)
      rc2_in = real(1/(gs/41.0), r_4)  ! transform µmol m-2 s-1 to s m-1
   end function canopy_resistence

   !=================================================================
   !=================================================================


   function stomatal_conductance(vpd_in,f1_in,g1) result(gs)
    ! return stomatal resistence based on Medlyn et al. 2011a
    ! Coded by Helena Alves do Prado

    use types, only: r_4 ,r_8
    use global_par, only: ca


    !implicit none

    real(r_4),intent(in) :: f1_in    !Photosynthesis (molCO2/m2/s)
    real(r_4),intent(in) :: vpd_in   !hPa
    real(r_4),intent(in) :: g1       ! model m (slope) (sqrt(kPa))
    real(r_8) :: gs       !Canopy conductance (molCO2 m-2 s-1)
    !     Internal
    !     --------
    real(r_8) :: D1       !sqrt(kPA)
    real(r_4) :: vapour_p_d

    vapour_p_d = vpd_in
    ! Assertions
    if(vpd_in .le. 0.0) vapour_p_d = 0.001
    if(vpd_in .gt. 4.0) vapour_p_d = 4.0
    ! print *, 'vpd going mad in canopy_resistence'
    ! stop
    ! endif

    D1 = sqrt(vapour_p_d)
    gs = 0.01 + 1.6 * (1.0 + (g1/D1)) * ((f1_in * 1e6)/ca)
    ! µmol m-2 s-1
 end function stomatal_conductance

 !=================================================================
 !=================================================================

   function water_ue(a, g, p0, vpd) result(wue)
      use types, only: r_4
      !implicit none

      real(r_4),intent(in) :: a, g, p0, vpd
      ! a = assimilacao; g = condutancia; p0 = pressao atm; vpd = vpd
      real(r_4) :: wue

      real(r_4) :: g_in, p0_in, e_in

      g_in = (1./g) * 41. ! convertendo a resistencia em condutancia mol m-2 s-1
      p0_in = p0 /10. ! convertendo pressao atm (mbar/hPa) em kPa
      e_in = g_in * (vpd/p0_in) ! calculando transpiracao

      if(a .eq. 0 .or. e_in .eq. 0) then
         wue = 0
      else
         wue = a/e_in
      endif
   end function water_ue

   !=================================================================
   !=================================================================

   function vapor_p_defcit(t,rh) result(vpd_0)
      use types, only: r_4
      !implicit none

      real(r_4),intent(in) :: t
      real(r_4),intent(in) :: rh

      real(r_4) :: es
      real(r_4) :: vpd_ac
      real(r_4) :: vpd_0

      ! ext func
      !real(r_4) :: tetens

      es = tetens(t)

      !     delta_e = es*(1. - ur)
      !VPD-REAL = Actual vapor pressure
      vpd_ac = es * rh       ! RESULT in hPa == mbar! we want kPa (DIVIDE by 10.)
      !Vapor Pressure Deficit
      vpd_0 = (es - vpd_ac) / 10.
   end function vapor_p_defcit

   !=================================================================
   !=================================================================

   function nrubisco(leaf_t,nbio_in) result(nb)
      use types
      use global_par, only: debug
      real(r_4), intent(in) :: leaf_t
      real(r_4), intent(in) :: nbio_in
      real(r_8) :: nb
      real(r_8) :: tl0, tlm, tl
      real(r_8) :: auxmax = 0.8, auxmin = 0.2

      tl0 = spec_leaf_area(1.0/12.0)
      tlm = spec_leaf_area(8.3)
      ! tlm --- 0.8
      ! tl0 --- 0.2
      if(leaf_t .lt. (1.0/12.0)) call abort
      if(leaf_t .gt. 8.3) call abort
      tl = spec_leaf_area(leaf_t)

      tl = 1.0-(((tl-tl0)/(tlm-tl0))*(auxmax-auxmin)+auxmin)

      if(debug) then
         if(tl < 0.0) print *, tl, 'total leaf Nrubisco - stop in nrubisco'
         if(tl < 0.0) stop
      endif
      nb = tl * nbio_in ! Leaf nitrogen that is rubisco

   end function nrubisco

   !=================================================================
   !=================================================================

   function nlignin(leaf_t, nbio_in) result(nl)
      use types
      use global_par, only: debug
      real(r_4), intent(in) :: leaf_t
      real(r_4), intent(in) :: nbio_in
      real(r_8) :: nl
      real(r_8) :: tl0, tlm, tl
      real(r_8) :: auxmax = 0.4, auxmin = 0.08

      tl0 = spec_leaf_area(1.0/12.0)
      tlm = spec_leaf_area(8.3)
      ! tlm --- 0.8
      ! tl0 --- 0.2
      tl = spec_leaf_area(leaf_t)

      tl =(((tl-tl0)/(tlm-tl0))*(auxmax-auxmin)+auxmin)

      if(debug) then
         if(tl < 0.0) print *, tl, 'total leaf Nlignin'
         if(tl < 0.0) stop
      endif
      nl = tl * nbio_in ! Leaf nitrogen that is rubisco

   end function nlignin

   !=================================================================
   !=================================================================

   subroutine photosynthesis_rate(temp,p0,ipar,ll,c4,nbio,pbio,&
        & leaf_turnover,storage_1,f1ab,vm,storage_2)

      ! f1ab SCALAR returns instantaneous photosynthesis rate at leaf level (molCO2/m2/s)
      ! vm SCALAR Returns maximum carboxilation Rate (Vcmax) (molCO2/m2/s)
      use types
      use utils
      use global_par
      use photo_par
      ! implicit none
      ! I
      real(r_4),intent(in) :: temp  ! temp °C
      real(r_4),intent(in) :: p0    ! atm Pressure hPa
      real(r_4),intent(in) :: ipar  ! mol Photons m-2 s-1
      real(r_4),intent(in) :: nbio  ! gm-2
      real(r_4),intent(in) :: pbio  ! gm-2
      logical(l_1),intent(in) :: ll ! is light limited?
      integer(i_4),intent(in) :: c4 ! is C4 Photosynthesis pathway?
      real(r_4),intent(in) :: leaf_turnover   ! y
      real(r_8),dimension(3),intent(in) :: storage_1 ! storage in vegetation (yesterday)
      ! O
      real(r_4),intent(out) :: f1ab ! Gross CO2 Assimilation Rate mol m-2 s-1
      real(r_8),intent(out) :: vm   ! PLS Vcmax mol m-2 s-1
      real(r_8),dimension(3),intent(out) :: storage_2 ! Storage in vegetation (now)

      real(r_8) :: f2,f3            !Michaelis-Menten CO2/O2 constant (Pa)
      real(r_8) :: mgama,vm_in      !Photo-respiration compensation point (Pa)
      real(r_8) :: rmax, r
      real(r_8) :: ci
      real(r_8) :: jp1
      real(r_8) :: jp2
      real(r_8) :: jp
      real(r_8) :: jc
      real(r_8) :: jl
      real(r_8) :: je,jcl
      real(r_8) :: b,c,c2,b2,es,j1,j2
      real(r_8) :: delta, delta2,aux_ipar
      real(r_8) :: f1a

      ! new vars C4 PHOTOSYNTHESIS
      real(r_8) :: ipar1
      real(r_8) :: tk           ! (K)
      real(r_8) :: t25          ! tk at 25°C (K)
      real(r_8) :: kp
      real(r_8) :: dummy0, dummy1, dummy2
      real(r_8) :: vpm, v4m
      real(r_8) :: cm, cm0, cm1, cm2

      real(r_8) :: vm_nutri
      real(r_8) :: nbio2, pbio2, xbio

      ! Calculating Fraction of leaf Nitrogen that is lignin
      xbio = nrubisco(leaf_turnover,nbio)
      !print *, xbio, "XBIO"
      ! n_lignin = nbio - xbio ! g m-2
      if(storage_1(2) .le. 0.0 ) then
         storage_2(2) = 0.0D0
         nbio2 = xbio
      else
         nbio2 = xbio + (storage_1(2) * 0.2)
         storage_2(2) = storage_1(2) - (storage_1(2) * 0.2)
      endif

      if(storage_1(3) .le. 0.0 ) then
         storage_2(3) = 0.0D0
         pbio2 = pbio
      else
         pbio2 = pbio + (storage_1(3) * 0.2)
         storage_2(3) = storage_1(3) - (storage_1(3) * 0.2)
      endif
      !print *, nbio2, 'nbio'
      !print *, pbio2, 'pbio'
      ! Saving storage values after a day of assimilation
      storage_2(1) = storage_1(1)

      !print*, "STORAGE1 "
      !print*, storage_1
      !print*, "Storage2"
      !print *, storage_2

      ! INCLUDING vcmax N and P Dependence
      ! Vmax dependence on N and P after Walker et al. 2014
      ! Prevent ilegal dlog application

      !print*, '--==--==--'
      !print*, nbio2, pbio2

      vm_nutri = 3.946D0 + (0.921D0 * dlog(nbio2)) - (0.121D0 * dlog(pbio2))
      vm_nutri = vm_nutri + (0.28D0 * dlog(nbio2) * dlog(pbio2))
      vm = dexp(vm_nutri) * 0.000001D0 ! Vcmax

      ! Vmax dependence on N and P after Domingues et al. 2010
      ! ... implement here.

      ! Rubisco Carboxilation Rate - temperature dependence
      vm_in = (vm*2.0**(0.1*(temp-25.0)))/(1.0+exp(0.3*(temp-36.0)))

      if(c4 .eq. 0) then
         !====================-C3 PHOTOSYNTHESIS-===============================
         !Photo-respiration compensation point (Pa)
         mgama = p3/(p8*(p9**(p10*(temp-p11))))
         !Michaelis-Menten CO2 constant (Pa)
         f2 = p12*(p13**(p10*(temp-p11)))
         !Michaelis-Menten O2 constant (Pa)
         f3 = p14*(p15**(p10*(temp-p11)))
         !Saturation vapour pressure (hPa)
         es = real(tetens(temp), r_8)
         !Saturated mixing ratio (kg/kg)
         rmax = 0.622*(es/(p0-es))
         !Moisture deficit at leaf level (kg/kg)
         r = -0.315*rmax
         !Internal leaf CO2 partial pressure (Pa)
         ci = p19* (1.-(r/p20)) * ((ca/9.901)-mgama) + mgama
         !Rubisco carboxilation limited photosynthesis rate (molCO2/m2/s)
         jc = vm_in*((ci-mgama)/(ci+(f2*(1.+(p3/f3)))))
         !Light limited photosynthesis rate (molCO2/m2/s)
         if (ll) then
            aux_ipar = ipar
         else
            aux_ipar = ipar - (ipar * 0.20)
         endif
         jl = p4*(1.0-p5)*aux_ipar*((ci-mgama)/(ci+(p6*mgama)))

         ! Transport limited photosynthesis rate (molCO2/m2/s) (RuBP) (re)generation
         ! ---------------------------------------------------
         je = p7*vm_in

         !Jp (minimum between jc and jl)
         !------------------------------
         b = (-1.)*(jc+jl)
         c = jc*jl
         delta = (b**2)-4.0*a*c
         ! if(delta .eq. 0.0)then
         !    jp = (-b) / (2 * a)
         ! else if(delta .gt. 0.0) then
         jp1 = (-b-(sqrt(delta)))/(2.0*a)
         jp2 = (-b+(sqrt(delta)))/(2.0*a)
         jp = amin1(jp1,jp2)
         ! else
         !    jp = 0.0
         ! endif

         !Leaf level gross photosynthesis (minimum between jc, jl and je)
         !---------------------------------------------------------------
         b2 = (-1.)*(jp+je)
         c2 = jp*je
         delta2 = (b2**2)-4.0*a2*c2
         ! if(delta2 .eq. 0.0)then
         !    f1a = (-b2) / (2.0 * a2)
         ! else if(delta2 .gt. 0.0) then
         j1 = (-b2-(sqrt(delta2)))/(2.0*a2)
         j2 = (-b2+(sqrt(delta2)))/(2.0*a2)
         f1a = amin1(j1,j2)
         ! else
         !    f1a = 0.0
         ! endif

         f1ab = real(f1a,r_4)
         return
      else
         !===========================-C4 PHOTOSYNTHESIS-=============================
         !  USE PHOTO_PAR
         ! ! from Chen et al. 1994
         tk = temp + 273.15           ! K
         t25 = 273.15 + 25.0          ! K
         kp = kp25 * (2.1**(0.1*(tk-t25))) ! ppm

         if (ll) then
            aux_ipar = ipar
         else
            aux_ipar = ipar - (ipar * 0.20)
         endif

         ipar1 = aux_ipar * 1e6  ! µmol m-2 s-1 - 1e6 converts mol to µmol

         !maximum PEPcarboxylase rate Arrhenius eq. (Dependence on temperature)
         dummy1 = 1.0 + exp((s_vpm * t25 - h_vpm)/(r_vpm * t25))
         dummy2 = 1.0 + exp((s_vpm * tk - h_vpm)/(r_vpm * tk))
         dummy0 = dummy1 / dummy2
         vpm =  vpm25 * exp((-e_vpm/r_vpm) * (1.0/tk - 1.0/t25)) * dummy0

         ! ! actual PEPcarboxylase rate under ipar conditions
         v4m = (alphap * ipar1) / sqrt(1 + alphap**2 * ipar1**2 / vpm**2)

         ! [CO2] mesophyl
         cm0 = 1.674 - 6.1294 * 10.0**(-2) * temp
         cm1 = 1.1688 * 10.0**(-3) * temp ** 2
         cm2 = 8.8741 * 10.0**(-6) * temp ** 3
         cm = 0.7 * ca * ((cm0 + cm1 - cm2) / 0.73)

         ! ! When light or PEP carboxylase is limiting
         ! ! FROM CHEN et al. 1994:
         jcl = ((V4m * cm) / (kp + cm)) * 1e-6   ! molCO2 m-2 s-1 / 1e-6 convets µmol 2 mol

         ! When V (RuBP regeneration) is limiting
         je = p7 * vm_in

         ! !Leaf level gross photosynthesis (minimum between jcl and je)
         ! !---------------------------------------------------------------
         b2 = (-1.)*(jcl+je)
         c2 = jcl*je
         delta2 = (b2**2)-4.0*a2*c2
         ! if(delta2 .eq. 0.0)then
         !    f1a = (-b2) / (2.0 * a2)
         ! else if(delta2 .gt. 0.0) then
         j1 = (-b2-(sqrt(delta2)))/(2.0*a2)
         j2 = (-b2+(sqrt(delta2)))/(2.0*a2)
         f1a = amin1(j1,j2)
         ! else
         !    f1a = 0.0
         ! endif

         f1ab = real(f1a,r_4)
         return
      endif
   end subroutine photosynthesis_rate

   !=================================================================
   !=================================================================

   !=====================================================================
   !c     subroutine allocation calculates the daily carbon content of each
   !c     compartment
   !c
   !c     code written by Bianca Rius & David Lapola (27.Ago.2015)
   !c     Modified 02-08-2017 jp - nutrient cycles implementation
   !c=====================================================================

   subroutine allocation(dt,npp,nmin,plab,scl1,sca1,scf1,storage,&
        &storage_out,scl2,sca2,scf2,leaf_litter,cwd,root_litter,&
        &nuptk,puptk,litter_nutrient_ratios,end_pls_day)

!!!!!&&&--- ALLOCATION
      use types
      use global_par, only: ntraits,npls,rfrac_leaf,rfrac_froot,rfrac_wood,cmin
      !implicit none

      ! variables I/O
      real(r_4),dimension(ntraits),intent(in) :: dt  ! PLS attributes
      real(r_4),intent(in) :: npp  ! npp (KgC/m2/yr) from assimilation process
      real(r_8),intent(in) :: scl1 ! previous day carbon content on leaf compartment (KgC/m2)
      real(r_8),intent(in) :: sca1 ! previous day carbon content on aboveground woody biomass compartment(KgC/m2)
      real(r_8),intent(in) :: scf1 ! previous day carbon content on fine roots compartment (KgC/m2)
      real(r_4),intent(in) :: nmin ! N in mineral N pool(g m-2)
      real(r_4),intent(in) :: plab ! P in labile pool (g m-2)
      real(r_8),dimension(3),intent(in) :: storage ! Three element array- storage pool([C,N,P]) g m-2
      ! O
      real(r_8),dimension(3),intent(out) :: storage_out
      real(r_8),intent(out) :: scl2 ! final carbon content on leaf compartment (KgC/m2)
      real(r_8),intent(out) :: sca2 ! final carbon content on aboveground woody biomass compartment (KgC/m2)
      real(r_8),intent(out) :: scf2 ! final carbon content on fine roots compartment (KgC/m2)
      real(r_8),intent(out) :: cwd  ! coarse wood debris (to litter)(C) g m-2
      real(r_8),intent(out) :: root_litter ! to litter g(C) m-2
      real(r_8),intent(out) :: leaf_litter ! to litter g(C) m-2
      real(r_8),intent(out) :: nuptk ! N plant uptake g(N) m-2
      real(r_8),intent(out) :: puptk ! P plant uptake g(P) m-2
      real(r_8),dimension(6),intent(out) :: litter_nutrient_ratios ! [(lln2c),(rln2c),(cwdn2c),(llp2c),(rlp2c),(cwdp2c)]
      logical(l_1), intent(out) :: end_pls_day ! ABORT MISSION SIGN
      ! internal variables
      real(r_8) :: scf2_tmp = 0.0D0 ! Store veg carbon pool in a 64bit fp
      real(r_8) :: sca2_tmp = 0.0D0
      real(r_8) :: scl2_tmp = 0.0D0
      real(r_8) :: npp_pot  = 0.0D0 ! potential npp g m-2 day-1
      real(r_8) :: n_npp_pot = 0.0D0, p_npp_pot = 0.0D0
      real(r_8) :: carbon_to_storage_p = 0.0D0, carbon_to_storage_n = 0.0D0

      real(r_8) :: npp_awood  = 0.0D0, npp_froot  = 0.0D0, npp_leaf  = 0.0D0 ! Partitioned npp (g(C) m-2 day-1)
      real(r_8) :: npp_awoodn = 0.0D0, npp_frootn = 0.0D0, npp_leafn = 0.0D0
      real(r_8) :: npp_awoodp = 0.0D0, npp_frootp = 0.0D0, npp_leafp = 0.0D0

      ! Auxiliary variables to calculate Plant Nutrient Uptake
      real(r_8) :: aux1 = 0.0D0, aux2 = 0.0D0, aux3 = 0.0D0
      real(r_8) :: nscl = 0.0D0  ! g(N) m-2
      real(r_8) :: nsca = 0.0D0  ! g(N) m-2
      real(r_8) :: nscf = 0.0D0  ! g(N) m-2
      real(r_8) :: pscl = 0.0D0  ! g(P) m-2
      real(r_8) :: psca = 0.0D0  ! g(P) m-2
      real(r_8) :: pscf = 0.0D0  ! g(P) m-2
      ! traits
      real(r_4) :: aleaf = 0.0     ! allocatation to plant compartments
      real(r_4) :: aawood = 0.0
      real(r_4) :: afroot = 0.0
      real(r_4) :: tleaf = 0.0     ! Residence time(yr)
      real(r_4) :: tawood = 0.0
      real(r_4) :: tfroot = 0.0
      real(r_4) :: leaf_n2c = 0.0  ! N:C ratios
      real(r_4) :: awood_n2c = 0.0
      real(r_4) :: froot_n2c = 0.0
      real(r_4) :: leaf_p2c = 0.0  ! P:C ratios
      real(r_4) :: awood_p2c = 0.0
      real(r_4) :: froot_p2c = 0.0

      real(r_8) :: avail_n = 0.0D0
      real(r_8) :: avail_p = 0.0D0
      real(r_8) :: potential_uptake_n = 0.0D0
      real(r_8) :: potential_uptake_p = 0.0D0

      ! Some flow control variables
      logical(l_1) :: no_cell = .false.
      logical(l_1) :: no_limit = .false.
      logical(l_1) :: n_limited = .false.
      logical(l_1) :: p_limited = .false.
      logical(l_1) :: no_allocation = .false.

      real(r_4) :: mult_factor = 1.0
      real(r_8) :: test_n_limitation
      real(r_8) :: test_p_limitation

      ! initialize some outputs
      ! If end_pls_day then PLS in not in the gridcell
      end_pls_day = .false.
      storage_out = (/0.0D0, 0.0D0, 0.0D0/)

      ! First check for the carbon content in leafs and fine roots (scl1 & scf1).
      ! A little C content in these pools means a rare strategy that remains in the system
      ! During the spinup we need to check if these pools are r(n)ea(r)ly in 'steady state'

      ! If there are no carbon in fine roots AND leafs:
      ! Then PLS 'Die' Label 10 to zero outputs and next pls. Means that the rest of this sub isn't executed
      if(((scf1 .lt. cmin) .and. (scl1 .lt. cmin))) then !
         no_cell = .true.
         end_pls_day = .true.
         goto 10
      endif

      ! Catch the functional/demographic traits of pls
      tleaf  = dt(3) ! RESIDENCE TIME years
      tawood = dt(4)
      tfroot = dt(5)
      aleaf  = dt(6) ! ALLOCATION
      aawood = dt(7)
      afroot = dt(8)

      leaf_n2c  = dt(10) ! Nutrient:Carbon Ratios
      awood_n2c = dt(11)
      froot_n2c = dt(12)
      leaf_p2c  = dt(13)
      awood_p2c = dt(14)
      froot_p2c = dt(15)

      ! Potential NPP.
      ! transform Kg m-2 yr-1 in g m-2 day-1 plus the C in storage pool.
      ! If there is not NPP then no allocation process! only deallocation label 294
      if(npp .le. 0.0 .and. storage(1) .le. 0.0) then
         no_allocation = .true.
         npp_awood = 0.0
         npp_froot = 0.0
         npp_leaf = 0.0
         nuptk = 0.0
         puptk = 0.0
         goto 294
      endif


      !# There is C to allocate
      no_cell = .false.
      no_allocation = .false.

      ! You have: kg m-2 year-1
      ! You want: g m-2 day-1
      npp_pot = (npp * 2.73791) ! Transform Kg m-2 Year-1 to g m-2 day
      if (npp_pot .lt. 0.0) npp_pot = 0.0
      if(storage(1) .gt. 0.0) then
         npp_pot = npp_pot + storage(1)
         storage_out(1) = 0.0D0
      endif
! According to the CAETÊ heuristic, NPP is partitioned into 3 fluxes (# 3 cveg pools). Each one goes to a
! specific pool. BUT the stoichiometry of each pls need to be satisfied i.e.
! N:C and P:C must remain the same in VEG carbon pools after the allocation process
! If there is not N and P in biodisponible pools then allocation of all
! carbon assimilated will be limited to mantain the pls N:P:C ratio. The carbon that is not
! allocated go to the storage pool. Otherwise, nutrient is taken from nmin and plab and
! all npp is allocated. nuptk and pout will catch the plant nutrient uptake realized by this PLS.
! Out of this function (in productivity.f90) nuptk and pout will be weighted by PLS abundance
! In the last part nutrient resorption is calculated. Resorbed nutrients go to the storage pool.
! Potential NPP for each compartment

      ! Partitioning NPP for CVEG pools
      npp_leaf = aleaf * npp_pot
      npp_froot = afroot * npp_pot
      npp_awood = aawood * npp_pot

      ! Nutrients needed in allocation process (Potential Plant nutrient uptake)
      nscl = npp_leaf * leaf_n2c    ! g(N) m-2
      nscf = npp_froot * froot_n2c  ! g(N) m-2
      nsca = npp_awood * awood_n2c  ! g(N) m-2

      pscl = npp_leaf * leaf_p2c    ! g(P) m-2
      pscf = npp_froot * froot_p2c  ! g(P) m-2
      psca = npp_awood * awood_p2c  ! g(P) m-2

      potential_uptake_n = nscl + nsca + nscf !g m⁻²
      potential_uptake_p = pscl + psca + pscf !g m⁻²

      ! Compare disponible nutrients with potential uptake
      ! Calculate available pools
      ! TODO - the influence of AM and Ptase are calculated before
      if(storage(2) .gt. 0.0) then
         avail_n = (mult_factor * nmin) + storage(2) !g m⁻²
         storage_out(2) = 0.0D0
      else
         avail_n = mult_factor * nmin
      endif

      if(storage(3) .gt. 0.0) then
         avail_p = (mult_factor * plab) + storage(3) !g m⁻²
         storage_out(3) = 0.0D0
      else
         avail_p = mult_factor * plab
      endif

      test_n_limitation = avail_n - potential_uptake_n ! g(N) m-2
      test_p_limitation = avail_p - potential_uptake_p ! g(P) m-2

      if(test_n_limitation .lt. 0.0) then
         n_limited = .true.
      else
         n_limited = .false.
      endif

      if(test_p_limitation .lt. 0.0)then
         p_limited = .true.
      else
         p_limited = .false.
      endif

      ! ! Then check for limitation
      ! if(.not. n_limited .and. .not. p_limited) then ! no LIMITATION
      !    ! Real nutrient Uptake
      !    nuptk = potential_uptake_n ! g(N) m-2
      !    puptk = potential_uptake_p ! g(P) m-2
      !    no_limit = .true.
      !    goto 294 ! GOTO deallocation
      ! endif

      ! Calculate N and P limitation-----------------------------------------------

      ! n/puptk is the quantity of N/P that cannot be allocated

      ! N limitation --------------------------------------------------------------
      ! TODO include the value of limitation as an output

      if(n_limited) then     ! If nuptk is negative
         ! Calculate the missing N for each CVEG pool
         aux1 = abs(test_n_limitation) * aleaf     ! g(N) m-2 - Nitrogen amount
         aux2 = abs(test_n_limitation) * aawood
         aux3 = abs(test_n_limitation) * afroot
      else
         ! There is no N limitation
         aux1 = 0.0D0
         aux2 = 0.0D0
         aux3 = 0.0D0
         n_npp_pot = npp_pot
         carbon_to_storage_n = 0.0D0
         goto 66
      endif

      ! CALCULATE THE AMOUNT OF C THAT GOES TO STORAGE (Assimilated but not allocated)
      aux1 = aux1 * (leaf_n2c**(-1))  ! g(C) m-2 in leaves that cannot be allocated
      aux2 = aux2 * (awood_n2c**(-1)) ! g(C) m-2 in awood that cannot be allocated
      aux3 = aux3 * (froot_n2c**(-1)) ! g(C) m-2 in froots that cannot be allocated

      !After calculation this C can be stored
      carbon_to_storage_n = aux1 + aux2 + aux3
      ! C to be allocated
      ! TOTAL NPP ALLOCATED GIVEN th Nitrogen Limitation
      n_npp_pot = npp_pot - carbon_to_storage_n

66    continue
      ! Calculating real npp
      npp_leafn = n_npp_pot * aleaf
      npp_awoodn = n_npp_pot * aawood
      npp_frootn = n_npp_pot * afroot

      ! END N limitation ------------------------------------------------------------

      ! P limitation ----------------------------------------------------------------
      ! Using the same approach of N limitation

      if(p_limited) then
         aux1 = abs(test_p_limitation) * aleaf           ! g(P) m-2 - Phosphorus limitation is weighted by
         aux2 = abs(test_p_limitation) * aawood          ! allocation coefficients
         aux3 = abs(test_p_limitation) * afroot
      else
         aux1 = 0.0
         aux2 = 0.0
         aux3 = 0.0
         p_npp_pot = npp_pot
         carbon_to_storage_p = 0.0D0
         goto 99
      endif

      ! CALCULATE THE AMOUNT OF C THAT GOES TO STORAGE (Assimilated but not allocated)
      aux1 = aux1 * (leaf_p2c**(-1))  ! g(C) m-2 in leaves that cannot be allocated
      aux2 = aux2 * (awood_p2c**(-1)) ! g(C) m-2 in awood that cannot be allocated
      aux3 = aux3 * (froot_p2c**(-1)) ! g(C) m-2 in froots that cannot be allocated

      carbon_to_storage_p = aux1 + aux2 + aux3
      ! C to be allocated
      ! TOTAL NPP ALLOCATED GIVEN th Nitrogen Limitation
      p_npp_pot = npp_pot - carbon_to_storage_p
99    continue

      ! real NPP
      npp_leafp =  p_npp_pot * aleaf
      npp_awoodp = p_npp_pot * aawood
      npp_frootp = p_npp_pot * afroot

      ! END P limitation-------------------------------------------------------------

      ! NPP = min(NPP_N, NPP_P) i.e min(p_npp_pot, n_npp_pot)
      ! If ALLOCATION is P limited
 !     if (n_limited .or. p_limited) then

         if(p_npp_pot .lt. n_npp_pot) then
            p_limited = .true.
            no_allocation = .false.
            no_limit = .false.
            no_cell = .false.
            !print*, "P LIMTED"
            npp_leaf = npp_leafp   ! g(C) m-2 day-1
            npp_awood = npp_awoodp ! g(C) m-2 day-1
            npp_froot = npp_frootp ! g(C) m-2 day-1

            puptk = avail_p ! g(P) m-2
            nuptk = (npp_leaf * leaf_n2c) + (npp_awood * awood_n2c) + (npp_froot * froot_n2c) ! g(N) m-2
            if(nuptk .gt. avail_n) nuptk = avail_n
            storage_out(1) = carbon_to_storage_p ! g(C) m-2

            ! Else ALLOCATION is N limited
         else if(n_npp_pot .lt. p_npp_pot) then
            n_limited = .true.
            no_allocation = .false.
            no_limit = .false.
            no_cell = .false.
            !print*, "N LIMTED"
            npp_leaf = npp_leafn   ! g(C) m-2 day-1
            npp_awood = npp_awoodn ! g(C) m-2 day-1
            npp_froot = npp_frootn ! g(C) m-2 day-1

            nuptk = avail_n ! g(N) m-2
            puptk = (npp_leaf * leaf_p2c) + (npp_awood * awood_p2c) + (npp_froot * froot_p2c) ! g(P) m-2
            if(puptk .gt. avail_p) puptk = avail_p
            storage_out(1) = carbon_to_storage_n ! g(C) m-2

         else
            no_allocation = .false.
            no_limit = .false.
            no_cell = .false.
            ! check for colimitation
            !print*, "COLIMITED"
            !if(abs(p_npp_pot - n_npp_pot) .lt. 0.00000001) then
               npp_leaf = (npp_leafn + npp_leafp) / 2.0D0    ! g(C) m-2 day-1
               npp_awood = (npp_awoodn + npp_awoodp) / 2.0D0 ! g(C) m-2 day-1
               npp_froot = (npp_frootn + npp_frootp) / 2.0D0 ! g(C) m-2 day-1
               nuptk = (npp_leaf * leaf_n2c) + (npp_awood * awood_n2c) + (npp_froot * froot_n2c) ! g(N) m-2
               puptk = (npp_leaf * leaf_p2c) + (npp_awood * awood_p2c) + (npp_froot * froot_p2c) ! g(P) m-2
               storage_out(1) = (carbon_to_storage_n + carbon_to_storage_p) /2.0D0 ! g(C) m-2
         endif
      ! else
      !    ! No Limitation
      !    ! Real nutrient Uptake
      !    nuptk = potential_uptake_n ! g(N) m-2
      !    puptk = potential_uptake_p ! g(P) m-2
      !    no_limit = .true.
      !     goto 294 ! GOTO deallocation
      ! endif


      ! DEALLOCATION PROCESS
      ! ## Make calculations only for leaf and froots that are commmom to all PLSs
      ! Carbon pool times Turnover Rate (inverse of residence time)
      ! Shit happens here (UNDERFLOW) .edit. NO MORE?

294   continue ! Material going to soil + updating veg pools

      if(no_allocation) then
         npp_awood = 0.0D0
         npp_froot = 0.0D0
         npp_leaf = 0.0D0
      endif

      if(no_limit) then
         storage_out(1) = 0.0D0
      endif

      root_litter = ((1D3 * scf1) * (tfroot * 365.242)**(-1)) !/ tfroot! g(C) m-2
      leaf_litter = ((1D3 * scl1) * (tleaf * 365.242)**(-1)) !/ tleaf ! g(C) m-2

      ! calculate the C content of each compartment
      scl2_tmp = (1D3 * scl1) + npp_leaf - leaf_litter
      scf2_tmp = (1D3 * scf1) + npp_froot - root_litter

      scf2 = real(scf2_tmp,r_4) ! g(C) m-2 day-1
      scl2 = real(scl2_tmp,r_4) ! g(C) m-2 day-1
      ! ## if it's a woody strategy:
      if(tawood .gt. 0.0 .and. aawood .gt. 0.0 .and. sca1 .gt. 0.0) then
         cwd = ((1D3 * sca1) * (tawood * 365.242)**(-1)) !/ tawood! g(C) m-2
         sca2_tmp = (1D3 * sca1) + npp_awood - cwd  ! g(C) m-2
         sca2 = real(sca2_tmp,r_4) ! results in g(C) m-2 day-1
      else
         sca2 = 0.0
         cwd = 0.0
      endif

      ! Nutrient resorption
      ! N resorbed
      aux1 = (leaf_litter * leaf_n2c) * rfrac_leaf   ! g(N) m-2
      aux2 = (root_litter * froot_n2c) * rfrac_froot ! g(N) m-2
      if(aawood .gt. 0.0) then
         aux3 = (cwd * awood_n2c) * rfrac_wood       ! g(N) m-2
      else
         aux3 = 0.0
      endif

      ! N resorbed goes to storage pool
      storage_out(2) = aux1 + aux2 + aux3         ! g(N) m-2

      aux1 =  (leaf_litter * leaf_n2c) - aux1    !nitrogen in litter
      aux2 =  (root_litter * froot_n2c) - aux2
      if(aawood .gt. 0.0) then
         aux3 = (cwd * awood_n2c) - aux3
      else
         aux3 = 0.0
      endif

      if(leaf_litter .gt. 0.0) then
         litter_nutrient_ratios(1) = aux1 / leaf_litter ! N:C litter ratio g(N) g(C)-1
      else
         litter_nutrient_ratios(1) = 0.0
      endif

      if(root_litter .gt. 0.0) then
         litter_nutrient_ratios(2) = aux2 / root_litter ! N:C litter ratio g(N) g(C)-1
      else
         litter_nutrient_ratios(2) = 0.0
      endif

      ! NOT EVERY PLS HAVE A WOOD POOL
      ! If it is a woody strategy
      if(aawood .gt. 0.0) then
         if(cwd .gt. 0.0) then
            litter_nutrient_ratios(3) = aux3 / cwd      ! N:C litter ratio g(N) g(C)-1
         else
             litter_nutrient_ratios(3) = 0.0D0
         endif
      else
         litter_nutrient_ratios(3) = 0.0D0
      endif


      ! P resobed
      aux1 = (leaf_litter * leaf_p2c) * rfrac_leaf
      aux2 = (root_litter * froot_p2c) * rfrac_froot
      if(aawood .gt. 0.0) then
         aux3 = (cwd * awood_p2c) * rfrac_wood
      else
         aux3 = 0.0
      endif

      storage_out(3) = aux1 + aux2 + aux3

      aux1 =  (leaf_litter * leaf_p2c) - aux1    !nitrogen in litter
      aux2 =  (root_litter * froot_p2c) - aux2
      if(aawood .gt. 0.0) then
         aux3 = (cwd * awood_p2c) - aux3
      else
         aux3 = 0.0
      endif

      if(leaf_litter .gt. 0.0) then
         litter_nutrient_ratios(4) = aux1 / leaf_litter ! P:C litter ratio g(P) g(C)-1
      else
         litter_nutrient_ratios(4) = 0.0
      endif

      if(root_litter .gt. 0.0) then
         litter_nutrient_ratios(5) = aux2 / root_litter ! P:C litter ratio g(P) g(C)-1
      else
         litter_nutrient_ratios(5) = 0.0
      endif
! NOT EVERY PLS HAVE A WOOD POOL
      ! If it is a woody strategy
      if(aawood .gt. 0.0) then
         if(cwd .gt. 0.0) then
            litter_nutrient_ratios(6) = aux3 / cwd      ! N:C litter ratio g(N) g(C)-1
         else
             litter_nutrient_ratios(6) = 0.0D0
         endif
      else
         litter_nutrient_ratios(6) = 0.0D0
      endif


      ! END CALCULATIONS

      no_cell = .false.
      end_pls_day = .false.

      if(.not. no_cell) then
         scl2 = scl2 * 1e-3
         scf2 = scf2 * 1e-3
         sca2 = sca2 * 1e-3
      endif

      if(no_allocation .and. .not. no_cell) then
         nuptk = 0.0
         puptk = 0.0
         end_pls_day = .false.
      endif

10    continue
      if(no_cell) then
         end_pls_day = .true.
         storage_out = 0.0
         scl2 = 0.0
         scf2 = 0.0
         sca2 = 0.0
         leaf_litter = 0.0
         cwd = 0.0
         root_litter = 0.0
         nuptk = 0.0
         puptk = 0.0
         litter_nutrient_ratios = 0.0
      endif
      if(nuptk .lt. 0.0) nuptk = 0.0
      if(puptk .lt. 0.0) puptk = 0.0

      return

   end subroutine allocation

   ! ===========================================================
   ! ===========================================================

   subroutine spinup3(nppot,dt,cleafini,cfrootini,cawoodini)
      use types
      implicit none

      !parameters
      integer(kind=i_4),parameter :: ntl=36525

      ! inputs
      integer(kind=i_4) :: kk, k

      real(kind=r_4),intent(in) :: nppot
      real(kind=r_4),dimension(6),intent(in) :: dt
      ! intenal
      real(kind=r_4) :: sensitivity
      real(kind=r_4) :: nppot2
      ! outputs
      real(kind=r_4),intent(out) :: cleafini
      real(kind=r_4),intent(out) :: cawoodini
      real(kind=r_4),intent(out) :: cfrootini

      ! more internal
      real(kind=r_4),dimension(ntl) :: cleafi_aux
      real(kind=r_4),dimension(ntl) :: cfrooti_aux
      real(kind=r_4),dimension(ntl) :: cawoodi_aux

      real(kind=r_4) :: aux_leaf
      real(kind=r_4) :: aux_wood
      real(kind=r_4) :: aux_root
      real(kind=r_4) :: out_leaf
      real(kind=r_4) :: out_wood
      real(kind=r_4) :: out_root

      real(kind=r_4) :: aleaf  !npp percentage alocated to leaf compartment
      real(kind=r_4) :: aawood !npp percentage alocated to aboveground woody biomass compartment
      real(kind=r_4) :: afroot !npp percentage alocated to fine roots compartmentc
      real(kind=r_4) :: tleaf  !turnover time of the leaf compartment (yr)
      real(kind=r_4) :: tawood !turnover time of the aboveground woody biomass compartment (yr)
      real(kind=r_4) :: tfroot !turnover time of the fine roots compartment
      logical(kind=l_1) :: iswoody

      ! catch 'C turnover' traits
      tleaf  = dt(1)
      tawood = dt(2)
      tfroot = dt(3)
      aleaf  = dt(4)
      aawood = dt(5)
      afroot = dt(6)

      iswoody = ((aawood .gt. 0.0) .and. (tawood .gt. 0.0))

      sensitivity = 1.001
      if(nppot .le. 0.0) goto 200
      nppot2 = nppot !/real(npls,kind=r_4)
      do k=1,ntl
         if (k.eq.1) then
            cleafi_aux (k) =  aleaf * nppot2
            cawoodi_aux(k) = aawood * nppot2
            cfrooti_aux(k) = afroot * nppot2
         else
            aux_leaf = cleafi_aux(k-1) + (aleaf * nppot2)
            aux_wood = cawoodi_aux(k-1) + (aleaf * nppot2)
            aux_root = cfrooti_aux(k-1) + (afroot * nppot2)

            out_leaf = aux_leaf - (cleafi_aux(k-1) / tleaf)
            out_wood = aux_wood - (cawoodi_aux(k-1) / tawood)
            out_root = aux_root - (cfrooti_aux(k-1) / tfroot)

            if(iswoody) then
               cleafi_aux(k) = amax1(0.0, out_leaf)
               cawoodi_aux(k) = amax1(0.0, out_wood)
               cfrooti_aux(k) = amax1(0.0, out_root)
            else
               cleafi_aux(k) = amax1(0.0, out_leaf)
               cfrooti_aux(k) = amax1(0.0, out_root)
               cawoodi_aux(k) = 0.0
            endif

            kk =  floor(k*0.66)
            if(iswoody) then
               if((cfrooti_aux(k)/cfrooti_aux(kk).lt.sensitivity).and.&
                    &(cleafi_aux(k)/cleafi_aux(kk).lt.sensitivity).and.&
                    &(cawoodi_aux(k)/cawoodi_aux(kk).lt.sensitivity)) then

                  cleafini = cleafi_aux(k) ! carbon content (kg m-2)
                  cfrootini = cfrooti_aux(k)
                  cawoodini = cawoodi_aux(k)
                  !  print *, 'woody exitet in', k
                  exit
               endif
            else
               if((cfrooti_aux(k)&
                    &/cfrooti_aux(kk).lt.sensitivity).and.&
                    &(cleafi_aux(k)/cleafi_aux(kk).lt.sensitivity)) then

                  cleafini = cleafi_aux(k) ! carbon content (kg m-2)
                  cfrootini = cfrooti_aux(k)
                  cawoodini = 0.0
                  !  print *, 'grass exitet in', k
                  exit
               endif
            endif
         endif
      enddo                  !nt
200   continue
   end subroutine spinup3

   ! ===========================================================
   ! ===========================================================

   subroutine spinup2(nppot,dt,cleafini,cfrootini,cawoodini)
      use types
      use global_par, only: ntraits,npls
      implicit none

      !parameters
      integer(kind=i_4),parameter :: ntl=36525

      ! inputs
      integer(kind=i_4) :: i6, kk, k

      real(kind=r_4),intent(in) :: nppot
      real(kind=r_4),dimension(ntraits, npls),intent(in) :: dt
      ! intenal
      real(kind=r_4) :: sensitivity
      real(kind=r_4) :: nppot2
      ! outputs
      real(kind=r_4),dimension(npls),intent(out) :: cleafini
      real(kind=r_4),dimension(npls),intent(out) :: cawoodini
      real(kind=r_4),dimension(npls),intent(out) :: cfrootini

      ! more internal
      real(kind=r_4),dimension(ntl) :: cleafi_aux
      real(kind=r_4),dimension(ntl) :: cfrooti_aux
      real(kind=r_4),dimension(ntl) :: cawoodi_aux

      real(kind=r_4) :: aux_leaf
      real(kind=r_4) :: aux_wood
      real(kind=r_4) :: aux_root
      real(kind=r_4) :: out_leaf
      real(kind=r_4) :: out_wood
      real(kind=r_4) :: out_root

      real(kind=r_4),dimension(npls) :: aleaf  !npp percentage alocated to leaf compartment
      real(kind=r_4),dimension(npls) :: aawood !npp percentage alocated to aboveground woody biomass compartment
      real(kind=r_4),dimension(npls) :: afroot !npp percentage alocated to fine roots compartmentc
      real(kind=r_4),dimension(npls) :: tleaf  !turnover time of the leaf compartment (yr)
      real(kind=r_4),dimension(npls) :: tawood !turnover time of the aboveground woody biomass compartment (yr)
      real(kind=r_4),dimension(npls) :: tfroot !turnover time of the fine roots compartment
      logical(kind=l_1) :: iswoody

      ! catch 'C turnover' traits
      tleaf  = dt(3,:)
      tawood = dt(4,:)
      tfroot = dt(5,:)
      aleaf  = dt(6,:)
      aawood = dt(7,:)
      afroot = dt(8,:)

      sensitivity = 1.01
      if(nppot .le. 0.0) goto 200
      nppot2 = nppot !/real(npls,kind=r_4)
      do i6=1,npls
         iswoody = ((aawood(i6) .gt. 0.0) .and. (tawood(i6) .gt. 0.0))
         do k=1,ntl
            if (k .eq. 1) then
               cleafi_aux (k) =  aleaf(i6) * nppot2
               cawoodi_aux(k) = aawood(i6) * nppot2
               cfrooti_aux(k) = afroot(i6) * nppot2

            else
               aux_leaf = cleafi_aux(k-1) + (aleaf(i6) * nppot2)
               aux_wood = cawoodi_aux(k-1) + (aleaf(i6) * nppot2)
               aux_root = cfrooti_aux(k-1) + (afroot(i6) * nppot2)

               out_leaf = aux_leaf - (cleafi_aux(k-1) / tleaf(i6))
               out_wood = aux_wood - (cawoodi_aux(k-1) / tawood(i6))
               out_root = aux_root - (cfrooti_aux(k-1) / tfroot(i6))

               if(iswoody) then
                  cleafi_aux(k) = amax1(0.0, out_leaf)
                  cawoodi_aux(k) = amax1(0.0, out_wood)
                  cfrooti_aux(k) = amax1(0.0, out_root)
               else
                  cleafi_aux(k) = amax1(0.0, out_leaf)
                  cawoodi_aux(k) = 0.0
                  cfrooti_aux(k) = amax1(0.0, out_root)
               endif

               kk =  floor(k*0.66)
               if(iswoody) then
                  if((cfrooti_aux(k)/cfrooti_aux(kk).lt.sensitivity).and.&
                       &(cleafi_aux(k)/cleafi_aux(kk).lt.sensitivity).and.&
                       &(cawoodi_aux(k)/cawoodi_aux(kk).lt.sensitivity)) then

                     cleafini(i6) = cleafi_aux(k) ! carbon content (kg m-2)
                     cfrootini(i6) = cfrooti_aux(k)
                     cawoodini(i6) = cawoodi_aux(k)
                     exit
                  endif
               else
                  if((cfrooti_aux(k)&
                       &/cfrooti_aux(kk).lt.sensitivity).and.&
                       &(cleafi_aux(k)/cleafi_aux(kk).lt.sensitivity)) then

                     cleafini(i6) = cleafi_aux(k) ! carbon content (kg m-2)
                     cfrootini(i6) = cfrooti_aux(k)
                     cawoodini(i6) = 0.0
                     exit
                  endif
               endif
            endif
         enddo                  !nt
      enddo                     !npls
200   continue
   end subroutine spinup2

  !===================================================================
  !===================================================================
!   CORRIGIR A TEPERATURA DO SOLO

   function m_resp(temp,cl1_mr,cf1_mr,ca1_mr,&
        & sto_mr,n2cl,n2cw,n2cf,aawood_mr) result(rm)
      use types, only: r_4,r_8
      use global_par, only: sapwood
      !implicit none

      real(r_4), intent(in) :: temp
      real(r_8), intent(in) :: cl1_mr
      real(r_8), intent(in) :: cf1_mr
      real(r_8), intent(in) :: ca1_mr
      real(r_8), dimension(3), intent(in) :: sto_mr
      real(r_4), intent(in) :: n2cl
      real(r_4), intent(in) :: n2cw
      real(r_4), intent(in) :: n2cf
      real(r_4), intent(in) :: aawood_mr
      real(r_4) :: rm

      real(r_8) :: csa, rm64, rml64,stoc,ston
      real(r_8) :: rmf64, rms64, storage_resp

      !   Autothrophic respiration
      !   ========================
      !   Maintenance respiration (kgC/m2/yr) (based in Ryan 1991)

      stoc = sto_mr(1)
      ston = sto_mr(2)

      if(stoc < 0.0) stoc = 0.0
      if(ston < 0.0) ston = 0.0

      ! sapwood carbon content (kgC/m2). X% of woody tissues (Pavlick, 2013)
      ! only for woody PLSs
      if(aawood_mr .gt. 0.0) then
         csa = sapwood * ca1_mr
         rms64 = ((n2cw * (csa * 1e3)) * 15.0 * exp(0.03*temp))
      else
         rms64 = 0.0
      endif

      rml64 = ((n2cl * (cl1_mr * 1e3)) * 15.0 * exp(0.03*temp))

      rmf64 = ((n2cf * (cf1_mr * 1e3)) * 15.0 * exp(0.03*temp))

      storage_resp = ((ston * stoc) * 15.0 * exp(0.03*temp))

      rm64 = (rml64 + rmf64 + rms64 + storage_resp) * 1e-3

      rm = real(rm64,r_4)

      if (rm .lt. 0) then
         rm = 0.0
      endif

   end function m_resp

   !====================================================================
   !====================================================================

   function g_resp(beta_leaf,beta_awood, beta_froot,aawood_rg) result(rg)
      use types, only: r_4,r_8
      !implicit none

      real(r_8), intent(in) :: beta_leaf
      real(r_8), intent(in) :: beta_froot
      real(r_8), intent(in) :: beta_awood
      real(r_4), intent(in) :: aawood_rg
      real(r_4) :: rg

      real(r_8) :: rg64, rgl64, rgf64, rgs64
      real(r_8) :: a1,a2,a3

      !     Autothrophic respiration
      !     Growth respiration (KgC/m2/yr)(based in Ryan 1991; Sitch et al.
      !     2003; Levis et al. 2004)

      a1 = beta_leaf
      a2 = beta_froot
      a3 = beta_awood

      if(a1 .le. 0.0) a1 = 0.0
      if(a2 .le. 0.0) a2 = 0.0
      if(a3 .le. 0.0) a3 = 0.0

      rgl64 = 1.25 * a1
      rgf64 = 1.25 * a2

      if(aawood_rg .gt. 0.0) then
         rgs64 = 1.25 * a3
      else
         rgs64 = 0.0
      endif

      rg64 = rgl64 + rgf64 + rgs64

      rg = real(rg64,r_4)

      if (rg.lt.0) then
         rg = 0.0
      endif

   end function g_resp

   !====================================================================
   !====================================================================

   function tetens(t) result(es)
      ! returns Saturation Vapor Pressure (hPa), using Buck equation

      ! buck equation...references:
      ! http://www.hygrometers.com/wp-content/uploads/CR-1A-users-manual-2009-12.pdf
      ! Hartmann 1994 - Global Physical Climatology p.351
      ! https://en.wikipedia.org/wiki/Arden_Buck_equation#CITEREFBuck1996

      ! Buck AL (1981) New Equations for Computing Vapor Pressure and Enhancement Factor.
      !      J. Appl. Meteorol. 20:1527–1532.

      use types, only: r_4
      !implicit none

      real(r_4),intent( in) :: t
      real(r_4) :: es

      if (t .ge. 0.) then
         es = 6.1121 * exp((18.729-(t/227.5))*(t/(257.87+t))) ! Arden Buck
         !es = es * 10 ! transform kPa in mbar == hPa
         return
      else
         es = 6.1115 * exp((23.036-(t/333.7))*(t/(279.82+t))) ! Arden Buck
         !es = es * 10 ! mbar == hPa ! mbar == hPa
         return
      endif

   end function tetens

   !====================================================================
   !====================================================================

   SUBROUTINE PFT_AREA_FRAC(CLEAF1, CFROOT1, CAWOOD1, OCP_COEFFS, OCP_WOOD)
      USE TYPES, ONLY: L_1, I_4, R_4, R_8
      USE GLOBAL_PAR, ONLY: NPLS, SAPWOOD
      !IMPLICIT NONE

      INTEGER(KIND=I_4),PARAMETER :: NPFT = NPLS ! PLSS FUTURAMENTE SERAO

      REAL(KIND=R_8),DIMENSION(NPFT),INTENT( IN) :: CLEAF1, CFROOT1, CAWOOD1
      REAL(KIND=R_8),DIMENSION(NPFT),INTENT(OUT) :: OCP_COEFFS
      LOGICAL(KIND=L_1),DIMENSION(NPFT),INTENT(OUT) :: OCP_WOOD

      REAL(KIND=R_8),DIMENSION(NPFT) :: CLEAF, CAWOOD, CFROOT
      REAL(KIND=R_8),DIMENSION(NPFT) :: TOTAL_BIOMASS_PFT,TOTAL_W_PFT
      INTEGER(KIND=I_4) :: P,I
      INTEGER(KIND=I_4),DIMENSION(1) :: MAX_INDEX
      REAL(KIND=R_8) :: TOTAL_BIOMASS, TOTAL_WOOD
      INTEGER(KIND=I_4) :: FIVE_PERCENT

      TOTAL_BIOMASS = 0.0
      TOTAL_WOOD = 0.0

      CLEAF = CLEAF1
      CFROOT = CFROOT1
      CAWOOD = CAWOOD1

      DO P = 1,NPFT
         TOTAL_W_PFT(P) = 0.0
         TOTAL_BIOMASS_PFT(P) = 0.0
         OCP_COEFFS(P) = 0.0
         OCP_WOOD(P) = .FALSE.
      ENDDO

      ! CHECK FOR NAN IN CLEAF CAWOOD CFROOT
      DO P = 1,NPFT
         IF(ISNAN(CLEAF(P))) CLEAF(P) = 0.0
         IF(ISNAN(CFROOT(P))) CFROOT(P) = 0.0
         IF(ISNAN(CAWOOD(P))) CAWOOD(P) = 0.0
      ENDDO


      DO P = 1,NPFT
         TOTAL_BIOMASS_PFT(P) = CLEAF(P) + CFROOT(P) + CAWOOD(P) ! (SAPWOOD * CAWOOD(P)) ! ONLY SAPWOOD?
         TOTAL_BIOMASS = TOTAL_BIOMASS + TOTAL_BIOMASS_PFT(P)
         TOTAL_WOOD = TOTAL_WOOD + CAWOOD(P)
         TOTAL_W_PFT(P) = CAWOOD(P)
      ENDDO

      !     GRID CELL OCCUPATION COEFFICIENTS
      IF(TOTAL_BIOMASS .GT. 0.0) THEN
         DO P = 1,NPFT
            OCP_COEFFS(P) = TOTAL_BIOMASS_PFT(P) / TOTAL_BIOMASS
            IF(OCP_COEFFS(P) .LT. 0.0) OCP_COEFFS(P) = 0.0
            !IF(ISNAN(OCP_COEFFS(P))) OCP_COEFFS(P) = 0.0
         ENDDO
      ELSE
         DO P = 1,NPFT
            OCP_COEFFS(P) = 0.0
         ENDDO
      ENDIF

      !     GRIDCELL PFT LIGTH LIMITATION BY WOOD CONTENT
      FIVE_PERCENT = NINT(REAL(NPFT) * 0.05)
      IF(FIVE_PERCENT .EQ. 0) FIVE_PERCENT = 1
      IF(FIVE_PERCENT .EQ. 1) THEN
         IF(TOTAL_WOOD .GT. 0.0) THEN
            MAX_INDEX = MAXLOC(TOTAL_W_PFT)
            I = MAX_INDEX(1)
            OCP_WOOD(I) = .TRUE.
         ENDIF
      ELSE
         DO P = 1,FIVE_PERCENT
            IF(TOTAL_WOOD .GT. 0.0) THEN
               MAX_INDEX = MAXLOC(TOTAL_W_PFT)
               I = MAX_INDEX(1)
               TOTAL_W_PFT(I) = 0.0
               OCP_WOOD(I) = .TRUE.
            ENDIF
         ENDDO
      ENDIF

   END SUBROUTINE PFT_AREA_FRAC

   !====================================================================
   !====================================================================

end module photo


module water

  ! this module defines functions related to surface water balance
  implicit none
  private

  ! functions defined here:

  public ::              &
       soil_temp        ,&
       soil_temp_sub    ,&
       penman           ,&
       evpot2           ,&
       available_energy ,&
       runoff


contains

  !=================================================================
  !=================================================================

  subroutine soil_temp_sub(temp, tsoil)
  ! Calcula a temperatura do solo. Aqui vamos mudar no futuro!
  ! a tsoil deve ter relacao com a et realizada...
  ! a profundidade do solo (H) e o coef de difusao (DIFFU) devem ser
  ! variaveis (MAPA DE SOLO?; agua no solo?)
  use types
  use global_par
  !implicit none
  integer(i_4),parameter :: m = 1095

  real(r_4),dimension(m), intent( in) :: temp ! future __ make temps an allocatable array
  real(r_4), intent(out) :: tsoil

  ! internal vars

  integer(i_4) :: n, k
  real(r_4) :: t0 = 0.0
  real(r_4) :: t1 = 0.0

  tsoil = -9999.0

  do n=1,m !run to attain equilibrium
     k = mod(n,12)
     if (k.eq.0) k = 12
     t1 = (t0*exp(-1.0/tau) + (1.0 - exp(-1.0/tau)))*temp(k)
     tsoil = (t0 + t1)/2.0
     t0 = t1
  enddo
  end subroutine soil_temp_sub

  !=================================================================
  !=================================================================

  function soil_temp(t0,temp) result(tsoil)
    use types
    use global_par, only: h, tau, diffu
    !implicit none

    real(r_4),intent( in) :: temp
    real(r_4),intent( in) :: t0
    real(r_4) :: tsoil

    real(r_4) :: t1 = 0.0

    t1 = (t0*exp(-1.0/tau) + (1.0 - exp(-1.0/tau)))*temp
    tsoil = (t0 + t1)/2.0
  end function soil_temp

  !=================================================================
  !=================================================================

  function penman (spre,temp,ur,rn,rc2) result(evap)
    use types, only: r_4
    use global_par, only: rcmin, rcmax
    use photo, only: tetens
    !implicit none


    real(r_4),intent(in) :: spre                 !Surface pressure (mbar)
    real(r_4),intent(in) :: temp                 !Temperature (°C)
    real(r_4),intent(in) :: ur                   !Relative humidity (0-1)
    real(r_4),intent(in) :: rn                   !Radiation balance (W/m2)
    real(r_4),intent(in) :: rc2                  !Canopy resistence (s/m)

    real(r_4) :: evap                            !Evapotranspiration (mm/day)
    !     Parameters
    !     ----------
    real(r_4) :: ra, h5, t1, t2, es, es1, es2, delta_e, delta
    real(r_4) :: gama, gama2


    ra = rcmin
    h5 = 0.0275               !mb-1

    !     Delta
    !     -----
    t1 = temp + 1.
    t2 = temp - 1.
    es1 = tetens(t1)       !Saturation partial pressure of water vapour at temperature T
    es2 = tetens(t2)

    delta = (es1-es2)/(t1-t2) !mbar/oC
    !
    !     Delta_e
    !     -------
    es = tetens (temp)
    delta_e = es*(1. - ur)    !mbar

    if ((delta_e.ge.(1./h5)-0.5).or.(rc2.ge.rcmax)) evap = 0.
    if ((delta_e.lt.(1./h5)-0.5).or.(rc2.lt.rcmax)) then
       !     Gama and gama2
       !     --------------
       gama  = spre*(1004.)/(2.45e6*0.622)
       gama2 = gama*(ra + rc2)/ra

       !     Real evapotranspiration
       !     -----------------------
       evap = (delta* rn + (1.20*1004./ra)*delta_e)/(delta+gama2) !W/m2
       evap = evap*(86400./2.45e6) !mm/day
       evap = amax1(evap,0.)  !Eliminates condensation
    endif
  end function penman

  !=================================================================
  !=================================================================

  function available_energy(temp) result(ae)
    use types, only: r_4
    !implicit none

    real(r_4),intent(in) :: temp
    real(r_4) :: ae

    ae = 2.895 * temp + 52.326 !from NCEP-NCAR Reanalysis data
  end function  available_energy

  !=================================================================
  !=================================================================

  function runoff(wa) result(roff)
    use types, only: r_4
    !implicit none

    real(r_4),intent(in) :: wa
    real(r_4):: roff

    !  roff = 38.*((w/wmax)**11.) ! [Eq. 10]
    roff = 11.5*((wa)**6.6) !from NCEP-NCAR Reanalysis data
  end function  runoff

  !=================================================================
  !=================================================================

  function evpot2 (spre,temp,ur,rn) result(evap)
    use types, only: r_4
    use global_par, only: rcmin, rcmax
    use photo, only: tetens
    !implicit none

    !Commments from CPTEC-PVM2 code
!    c Entradas
!c --------
!c spre   = pressao aa supeficie (mb)
!c temp   = temperatura (oC)
!c ur     = umidade relativa  (0-1,adimensional)
!c rn     = saldo de radiacao (W m-2)
!c
!c Saida
!c -----
!c evap  = evapotranspiracao potencial sem estresse (mm/dia)

    !     Inputs

    real(r_4),intent(in) :: spre                 !Surface pressure (mb)
    real(r_4),intent(in) :: temp                 !Temperature (oC)
    real(r_4),intent(in) :: ur                   !Relative humidity (0-1,dimensionless)
    real(r_4),intent(in) :: rn                   !Radiation balance (W/m2)
    !     Output
    !     ------
    !
    real(r_4) :: evap                 !Evapotranspiration (mm/day)
    !     Parameters
    !     ----------
    real(r_4) :: ra, t1, t2, es, es1, es2, delta_e, delta
    real(r_4) :: gama, gama2, rc

    ra = rcmin            !s/m

    !     Delta

    t1 = temp + 1.
    t2 = temp - 1.
    es1 = tetens(t1)
    es2 = tetens(t2)
    delta = (es1-es2)/(t1-t2) !mb/oC

    !     Delta_e
    !     -------

    es = tetens (temp)
    delta_e = es*(1. - ur)    !mb

    !     Stomatal Conductance
    !     --------------------

    rc = rcmin

    !     Gama and gama2
    !     --------------

    gama  = spre*(1004.)/(2.45e6*0.622)
    gama2 = gama*(ra + rc)/ra

    !     Potencial evapotranspiration (without stress)
    !     ---------------------------------------------

    evap =(delta*rn + (1.20*1004./ra)*delta_e)/(delta+gama2) !W/m2
    evap = evap*(86400./2.45e6) !mm/day
    evap = amax1(evap,0.)     !Eliminates condensation
  end function evpot2

  !=================================================================
  !=================================================================

end module water
