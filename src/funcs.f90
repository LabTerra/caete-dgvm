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

! contacts :: David Montenegro Lapola <lapoladm ( at ) gmail.com>

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
        xylem_conductance      ,& ! (f), hydraulic conductance of xylem
        water_stress_modifier  ,& ! (f), F5 - water stress modifier (dimensionless)
        leaf_age_factor        ,& ! (f), effect of leaf age on photosynthetic rate
        photosynthesis_rate    ,& ! (s), leaf level CO2 assimilation rate (molCO2 m-2 s-1)
        canopy_resistence      ,& ! (f), Canopy resistence (from Medlyn et al. 2011a) (s/m) == m s-1
        stomatal_conductance   ,&
        ! scarbon_decaiment      ,& ! (f), Carbon decay borrowed from Pavlick et al. 2012
        vapor_p_defcit         ,& ! (f), Vapor pressure defcit  (kPa)
        tetens                 ,& ! (f), Maximum vapor pressure (hPa)
        nrubisco               ,& ! (f), Fraction of N not in lignin (disponible to rubisco)
        nlignin                ,& ! (f), Fraction of N in lignin
        m_resp                 ,& ! (f), maintenance respiration (plants)
        spinup2                ,& ! (s)
        spinup3                ,& ! (s)
        g_resp                 ,& ! (f), growth Respiration (kg m-2 yr-1)
        ! carbon2                ,& ! (s)
        ! carbon3                ,& ! (s), soil + litter + heterothrophic respiration
        pft_area_frac          ,& ! (s), area fraction by biomass
        pft_par                ,& ! (s), aux subroutine to read pls data               (D)
        ascii2bin              ,& ! (s), converts txt file (traits,lines) to .bin file (D)
        allocation             ,& ! (s), daily npp allocation subroutine
        water_ue               ,& ! (f), water use efficiency
        rc_canopy_scaling      ,& ! (f), Try it!
        leap                      ! (f), is leap year ?

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

   function leaf_area_index(cleaf,sla) result(lai)
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

      sunlai = (1.0-(exp(-p26*lai)))/p26
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
         lai_ss = real((1.0-(exp(-p26*sunlai)))/p26,r_4) !sun decl 90 degrees
         return
      endif

      if(fs .eq. 2) then
         !f4shade
         lai_ss = real((1.0-(exp(-p27*shadelai)))/p27,r_4) !sun decl ~20 degrees
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

   function xylem_conductance(psi_soil, psi_g, p50, vulnerability_curve) result(v)    ! based in Eller et al. 2018
      use types
      use global_par

      real(r_4), intent(in) :: psi_soil, psi_g, p50, vulnerability_curve
      real(r_4) :: V

      v = 1.0 / (1.0 + ((psi_soil - psi_g) / p50) ** vulnerability_curve)

   end function

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
   end function water_stress_modifier

   ! =============================================================
   ! =============================================================

   function leaf_age_factor(mu, acrit, a) result(fa)    ! based in Caldararu et al. 2018
      use types

      real(r_4), intent(in) :: mu, acrit, a
      real(r_4) :: fa

      fa = amin1(1.0, exp(mu * (acrit-a)))

   end function

   !=================================================================
   !=================================================================

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
      ! n_lignin = nbio - xbio ! g m-2
      nbio2 = xbio + (storage_1(2) * 0.5) ! Leaf nitrogen that is 'rubisco'
      pbio2 = pbio + (storage_1(3) * 0.5)

      ! print *, nbio2, 'nbio'
      ! print *, pbio2, 'pbio'
      ! Saving storage values after a day of assimilation
      storage_2(1) = storage_1(1)
      storage_2(2) = storage_1(2) - (storage_1(2) * 0.5)
      storage_2(3) = storage_1(3) - (storage_1(3) * 0.5)

      ! INCLUDING vcmax N and P Dependence
      ! Vmax dependence on N and P after Walker et al. 2014
      vm_nutri = 3.946 + 0.921 * log(nbio2) - 0.121 * log(pbio2)
      vm_nutri = vm_nutri + 0.282 * log(nbio2) * log(pbio2)
      vm = exp(vm_nutri) * 1e-6 ! Vcmax

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
      real(r_8),dimension(3), intent(in) :: scl1 ! previous day carbon content on leaf compartment (KgC/m2)
      real(r_8),intent(in) :: sca1 ! previous day carbon content on aboveground woody biomass compartment(KgC/m2)
      real(r_8),intent(in) :: scf1 ! previous day carbon content on fine roots compartment (KgC/m2)
      real(r_8),intent(in) :: nmin ! N in mineral N pool(kg m-2)
      real(r_8),intent(in) :: plab ! P in labile pool (kg m-2)
      real(r_8),dimension(3),intent(in) :: storage ! Three element array- storage pool([C,N,P]) g m-2
      ! O
      real(r_8),dimension(3),intent(out) :: storage_out
      real(r_8),dimension(3),intent(out) :: scl2 ! final carbon content on leaf compartment (KgC/m2)
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
      real(r_8) :: scf2_128 = 0.0 ! Store veg carbon pool in a 64bit fp
      real(r_8) :: sca2_128 = 0.0
      real(r_8), dimension(3) :: scl2_128 = 0.0
      real(r_8) :: npp_pot  = 0.0 ! potential npp g m-2 day-1

      real(r_8) :: noutn, noutp ! plant N/P uptake given Nutrient Limitation (g(N\P) m-2)
      real(r_8) :: poutn, poutp
      real(r_8) :: npp_awood, npp_froot, npp_leaf ! Partitioned npp (g(C) m-2 day-1)
      real(r_8) :: npp_awoodn, npp_frootn, npp_leafn
      real(r_8) :: npp_awoodp, npp_frootp, npp_leafp
      real(r_8) :: total_n, total_p, leafscl2, leafscl3! total npp given N or P limitation g(C) m-2 day-1

      ! Auxiliary variables to calculate Plant Nutrient Uptake
      real(r_8) :: aux1, aux2, aux3
      real(r_8) :: nscl  ! g(N) m-2
      real(r_8) :: nsca  ! g(N) m-2
      real(r_8) :: nscf  ! g(N) m-2
      real(r_8) :: pscl  ! g(P) m-2
      real(r_8) :: psca  ! g(P) m-2
      real(r_8) :: pscf  ! g(P) m-2
      ! traits
      real(r_4) :: aleaf     ! allocatation to plant compartments
      real(r_4) :: aawood
      real(r_4) :: afroot
      real(r_4) :: tleaf     ! Residence time(yr)
      real(r_4) :: tawood
      real(r_4) :: tfroot
      real(r_4) :: leaf_n2c  ! N:C ratios
      real(r_4) :: awood_n2c
      real(r_4) :: froot_n2c
      real(r_4) :: leaf_p2c  ! P:C ratios
      real(r_4) :: awood_p2c
      real(r_4) :: froot_p2c

      ! Some flow control variables
      logical(l_1) :: no_cell = .false.
      logical(l_1) :: no_limit = .false.
      logical(l_1) :: n_limited = .false.
      logical(l_1) :: p_limited = .false.
      logical(l_1) :: no_allocation = .false.

      ! initialize last output
      end_pls_day = .false.

      ! First check for the carbon content in leafs and fine roots (scl1 & scf1).
      ! A little C content in these pools means a rare strategy that remains in the system
      ! During the spinup we need to check if these pools are r(n)ea(r)ly in 'steady state'

      ! If there are no carbon in fine roots AND leafs:
      ! Then PLS 'Die' Label 10 to zero outputs and next pls. Means that the rest of this sub isn't executed
      if(((scf1 .lt. cmin) .and. (sum(scl1) .lt. cmin))) then !
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
      if(npp .le. 0.000001 ) then!.and. storage(1) .le. 0.0) then
         no_allocation = .true.
         npp_awood = 0.0
         npp_froot = 0.0
         npp_leaf = 0.0
         goto 294
      endif

      npp_pot = (npp * 2.73791)
      if(storage(1) .ge. 0.0) npp_pot = (npp * 2.73791) + storage(1)
      no_cell = .false.
      no_allocation = .false.

      ! According the CAETÊ NPP is partitioned into 3 fluxes (# 3 cveg pools). Each one goes to a
      ! specific pool. BUT the stoichiometry of each pls need to be satisfied i.e.
      ! N:C and P:C must remain the same in VEG carbon pools after the allocation process
      ! If there is not N and P in biodisponible pools then allocation of all
      ! carbon assimilated will be limited to mantain the pls N:P:C ratio. The carbon that is not
      ! allocated go to the storage pool. Otherwise, nutrient is taken from nmin and plab and
      ! all npp is allocated. nuptk and pout will catch the plant nutrient uptake realized by this PLS.
      ! Out of this function (in productivity.f90) nuptk and pout will be weighted by PLS abundance
      ! In the last part nutrient resorption is calculated. Resorbed nutrients go to the storage pool.
      ! Potential NPP for each compartment
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

      ! Compare disponible nutrients with potential uptake
      nuptk = ((nmin * 1e3) + storage(2)) - (nscl + nsca + nscf) ! g(N) m-2
      puptk = ((plab * 1e3) + storage(3)) - (pscl + psca + pscf) ! g(P) m-2

      ! todo catch the infinities and NaNs
      if(nuptk .eq. -0.00000000) nuptk = 0.0
      if(puptk .eq. -0.00000000) puptk = 0.0
      if(nuptk .lt. 0.0) n_limited = .true.
      if(puptk .lt. 0.0) p_limited = .true.

      ! Then check for limitation

      if(.not. n_limited .and. .not. p_limited) then ! no LIMITATION
         ! Real nutrient Uptake
         nuptk = nscl + nsca + nscf ! g(N) m-2
         puptk = pscl + psca + pscf ! g(P) m-2
         no_limit = .true.
         ! ! ----------------------------
         !  print *, '   '
         !  print *, nout, ' nout  '
         !  print *, pout, ' pout  '
         !  print *, ' goint to 100   '
         ! !-----------------------------
         goto 294 ! GOTO deallocation process or ...
      endif

      ! Calculate N and P limitation-----------------------------------------------

      ! nout/pout is the quantity of N/P that cannot be allocated

      ! N limitation --------------------------------------------------------------
      if(n_limited) then     ! If nout is negative
         aux1 = nuptk * aleaf     ! g(N) m-2 - Nitrogen limitation is weighted by
         aux2 = nuptk * aawood    ! allocation coefficients. These numbers must be < 0
         aux3 = nuptk * afroot

         if (aux1 .eq. -0.0000000) aux1 = 0.0
         if (aux2 .eq. -0.0000000) aux2 = 0.0
         if (aux3 .eq. -0.0000000) aux3 = 0.0
         ! ----------------------------
         !  print *, '   '
         !  print *, aux1, ' nout leaf  '
         !  print *, aux2, ' nout awood '
         !  print *, aux3, ' nout froot '
         ! ! -----------------------------
      else
         aux1 = 0.0
         aux2 = 0.0
         aux3 = 0.0
         goto 66
      endif

      ! need to know the carbon amount that cannot be allocated given my N limitation
      ! auxi need to be zero or a negative number.
      ! Check if it's nan
      if(isnan(aux1)) aux1 = 0.0
      if(isnan(aux2)) aux2 = 0.0
      if(isnan(aux3)) aux3 = 0.0
      ! ! Check if it's inf
      if(aux1 .eq. aux1 - 1) aux1 = 0.0
      if(aux2 .eq. aux2 - 1) aux2 = 0.0
      if(aux3 .eq. aux3 - 1) aux3 = 0.0

      if(leaf_n2c .gt. 0.0) then
         aux1 = aux1 * leaf_n2c**(-1)  ! g(C) m-2 in leaves that cannot be allocated
      else
         aux1 = 0.0
      endif

      if(awood_n2c .gt. 0.0) then
         aux2 = aux2 * awood_n2c**(-1) ! g(C) m-2 in awood that cannot be allocated
      else
         aux2 = 0.0
      endif

      if(froot_n2c .gt. 0.0) then
         aux3 = aux3 * froot_n2c**(-1) ! g(C) m-2 in froots that cannot be allocated
      else
         aux3 = 0.0
      endif

      ! ! ----------------------------
      ! print *, '   '
      ! print *, aux1, ' nalloc leaf  '
      ! print *, aux2, ' nalloc awood '
      ! print *, aux3, ' nalloc froot '
      ! print *, '   '
      ! ! -----------------------------

      ! Check if it's nan
      if(isnan(aux1)) aux1 = 0.0
      if(isnan(aux2)) aux2 = 0.0
      if(isnan(aux3)) aux3 = 0.0
      ! ! Check if it's inf
      if(aux1 .eq. aux1 - 1) aux1 = 0.0
      if(aux2 .eq. aux2 - 1) aux2 = 0.0
      if(aux3 .eq. aux3 - 1) aux3 = 0.0

66    continue
      ! Calculating real npp
      npp_leafn = npp_leaf + aux1
      npp_awoodn = npp_awood + aux2
      npp_frootn = npp_froot + aux3

      if(isnan(npp_leafn)) npp_leafn = 0.0
      if(isnan(npp_awoodn)) npp_awoodn = 0.0
      if(isnan(npp_frootn)) npp_frootn = 0.0

      if(npp_leafn .eq. npp_leafn - 1) npp_leafn = 0.0
      if(npp_awoodn .eq. npp_awoodn - 1) npp_awoodn = 0.0
      if(npp_frootn .eq. npp_frootn - 1) npp_frootn = 0.0

      ! ! ----------------------------
      ! print *, '   '
      ! print *, npp_leafn, ' npp_leaf  '
      ! print *, npp_awoodn, ' npp_awood '
      ! print *, npp_frootn, ' npp_froot '
      ! print *, '   '
      ! ! -----------------------------

      total_n = npp_leafn + npp_awoodn + npp_frootn

      ! Real Nutrient uptake
      if(.not. p_limited) then
         poutn = pscl + psca + pscf
      else
         poutn = plab * 1e-3
      endif

      noutn = nmin * 1e-3

      ! print *, poutn, 'poutn'
      ! print *, noutn, 'noutn'

      ! END N limitation ------------------------------------------------------------

      ! P limitation ----------------------------------------------------------------
      ! Using the same approach of N limitation
      ! I need to know the carbon amount that cannot be alnoutplocated given my P limitation
      ! auxi need to be zero or a negative number.

      if(p_limited) then
         aux1 = puptk * aleaf           ! g(P) m-2 - Phosphorus limitation is weighted by
         aux2 = puptk * aawood          ! allocation coefficients
         aux3 = puptk * afroot

         if (aux1 .eq. -0.0000000) aux1 = 0.0
         if (aux2 .eq. -0.0000000) aux2 = 0.0
         if (aux3 .eq. -0.0000000) aux3 = 0.0
         ! ! ----------------------------
         ! !  print *, '   '
         !  print *, aux1, ' pout leaf  '
         !  print *, aux2, ' pout awood '
         !  print *, aux3, ' pout froot '
         ! ! -----------------------------
      else
         aux1 = 0.0
         aux2 = 0.0
         aux3 = 0.0
         goto 99
      endif

      ! Check if it's nan
      if(isnan(aux1)) aux1 = 0.0
      if(isnan(aux2)) aux2 = 0.0
      if(isnan(aux3)) aux3 = 0.0
      ! ! Check if it's inf
      if(aux1 .eq. aux1 - 1) aux1 = 0.0
      if(aux2 .eq. aux2 - 1) aux2 = 0.0
      if(aux3 .eq. aux3 - 1) aux3 = 0.0

      ! I need to know the carbon amount that cannot be allocated given my P limitation
      if(leaf_p2c .gt. 0.0) then
         aux1 = aux1 * leaf_p2c**(-1)  ! g(C) m-2 in leaves that cannot be allocated
      else
         aux1 = 0.0
      endif

      if(awood_p2c .gt. 0.0) then
         aux2 = aux2 * awood_p2c**(-1) ! g(C) m-2 in awood that cannot be allocated
      else
         aux2 = 0.0
      endif

      if(froot_p2c .gt. 0.0) then
         aux3 = aux3 * froot_p2c**(-1) ! g(C) m-2 in froots that cannot be allocated
      else
         aux3 = 0.0
      endif

      ! aux1 = aux1 * leaf_p2c**(-1)  ! g(C) m-2 in leaf that cannot be allocated
      ! aux2 = aux2 * awood_p2c**(-1) ! g(C) m-2 in awood
      ! aux3 = aux3 * froot_p2c**(-1) ! g(C) m-2 in froots

      ! ! ----------------------------
      ! print *, '   '
      ! print *, aux1, ' palloc leaf  '
      ! print *, aux2, ' palloc awood '
      ! print *, aux3, ' palloc froot '
      ! print *, '   '
      ! ! -----------------------------

      ! Check if it's nan
      if(isnan(aux1)) aux1 = 0.0
      if(isnan(aux2)) aux2 = 0.0
      if(isnan(aux3)) aux3 = 0.0
      ! ! Check if it's inf
      if(aux1 .eq. aux1 - 1) aux1 = 0.0
      if(aux2 .eq. aux2 - 1) aux2 = 0.0
      if(aux3 .eq. aux3 - 1) aux3 = 0.0
99    continue
      ! real NPP
      npp_leafp = npp_leaf + aux1
      npp_awoodp = npp_awood + aux2
      npp_frootp = npp_froot + aux3

      if(isnan(npp_leafp)) npp_leafp = 0.0
      if(isnan(npp_awoodp)) npp_awoodp = 0.0
      if(isnan(npp_frootp)) npp_frootp = 0.0

      if(npp_leafp .eq. npp_leafp - 1) npp_leafp= 0.0
      if(npp_awoodp .eq. npp_awoodp - 1) npp_awoodp = 0.0
      if(npp_frootp .eq. npp_frootp - 1) npp_frootp = 0.0

      total_p = npp_leafp + npp_awoodp + npp_frootp

      ! Real Nutrient uptake
      if(.not. n_limited) then
         noutp = nscl + nsca + nscf
      else
         noutp = nmin * 1e-3
      endif

      poutp = plab * 1e-3
      ! END P limitation-------------------------------------------------------------

      ! NPP = min(NPP_N, NPP_P) i.e min(total_p, total_n)
      ! If npp is P limited
      if(total_p .lt. total_n) then
         npp_leaf = npp_leafp   ! g(C) m-2 day-1
         npp_awood = npp_awoodp ! g(C) m-2 day-1
         npp_froot = npp_frootp ! g(C) m-2 day-1
         nuptk = noutp           ! g(N) m-2
         puptk = poutp            ! g(P) m-2
         storage_out(1) = npp_pot - total_p ! g(C) m-2
         no_limit = .false.
         ! Else npp is N limited
      else if(total_p .gt. total_n) then
         npp_leaf = npp_leafn   ! g(C) m-2 day-1
         npp_awood = npp_awoodn ! g(C) m-2 day-1
         npp_froot = npp_frootn ! g(C) m-2 day-1
         nuptk = noutn            ! g(N) m-2
         puptk = poutn           ! g(P) m-2
         storage_out(1) = npp_pot - total_n ! g(C) m-2
         no_limit = .false.
      else
         ! colimitation
         npp_leaf = (npp_leafn + npp_leafp) / 2.    ! g(C) m-2 day-1
         npp_awood = (npp_awoodn + npp_awoodp) / 2. ! g(C) m-2 day-1
         npp_froot = (npp_frootn + npp_frootp) / 2. ! g(C) m-2 day-1
         nuptk = nmin * 1e-3           ! g(N) m-2
         puptk = plab * 1e-3           ! g(P) m-2
         storage_out(1) = npp_pot - ((total_n + total_p)/2) ! g(C) m-2
         no_limit = .false.
      endif


      ! (DE)ALLOCATION PROCESS(ES)
      ! ## Make calculations only for leaf and froots that are commmom to all PLSs
      ! Carbon pool times Turnover Rate (inverse of residence time)
      ! Shit happens here (UNDERFLOW) .edit. NO MORE?

294   continue ! Material going to soil + updating veg pools
      if(no_allocation) then
         npp_awood = 0.0
         npp_froot = 0.0
         npp_leaf = 0.0
      endif

      if(no_limit) then
         storage_out(1) = 0.0
      endif

      leaf_litter = ((1e3 * scl1(3)) * (tleaf * 365.242)**(-1)) !/ tleaf ! g(C) m-2
      leafscl2 = ((1e3 * scl1(1)) * (tleaf * 365.242)**(-1))
      leafscl3 = ((1e3 * scl1(3)) * (tleaf * 365.242)**(-1))
      ! calculate the C content of each compartment
      scl2_128(1) = (1e3 * scl1(1)) + npp_leaf - leafscl2
      scl2_128(2) = (1e3 * scl1(2)) + leafscl2 - leafscl3
      scl2_128(3) = (1e3 * scl1(3)) - leaf_litter + leafscl3
      scf2_128 = (1e3 * scf1) + npp_froot - root_litter

      root_litter = ((1e3 * scf1) * (tfroot * 365.242)**(-1)) !/ tfroot! g(C) m-2
      scf2_128 = (1e3 * scf1) + npp_froot - root_litter

      scf2 = real(scf2_128,r_4) ! g(C) m-2 day-1
      scl2 = real(scl2_128,r_4) ! g(C) m-2 day-1
      ! ## if it's a woody strategy:
      if(tawood .gt. 0.0 .and. aawood .gt. 0.0 .and. sca1 .gt. 0.0) then
         cwd = ((1e3 * sca1) * (tawood * 365.242)**(-1)) !/ tawood! g(C) m-2
         sca2_128 = (1e3 * sca1) + npp_awood - cwd  ! g(C) m-2
         sca2 = real(sca2_128,r_4) ! results in g(C) m-2 day-1
      else
         sca2 = 0.0
         cwd = 0.0
      endif

      ! print *, root_litter, 'root litter'
      ! print *, leaf_litter, 'leaf_litter'
      ! print *, cwd, 'cwd'
      ! print *, scf2, 'scf2'
      ! print *, scl2, 'scl2'
      ! print *, sca2, 'sca2'

      ! Nutrient resorption

      ! N resorbed
      aux1 = (leaf_litter * leaf_n2c) * rfrac_leaf   ! g(N) m-2
      aux2 = (root_litter * froot_n2c) * rfrac_froot ! g(N) m-2
      if(aawood .gt. 0.0) then
         aux3 = (cwd * awood_n2c) * rfrac_wood       ! g(N) m-2
      else
         aux3 = 0.0
      endif
      ! Check if it's nan
      if(isnan(aux1)) aux1 = 0.0
      if(isnan(aux2)) aux2 = 0.0
      if(isnan(aux3)) aux3 = 0.0
      ! ! Check if it's inf
      if(aux1 .eq. aux1 - 1) aux1 = 0.0
      if(aux2 .eq. aux2 - 1) aux2 = 0.0
      if(aux3 .eq. aux3 - 1) aux3 = 0.0

      storage_out(2) = aux1 + aux2 + aux3         ! g(N) m-2

      aux1 =  (leaf_litter * leaf_n2c) - aux1    !nitrogen in litter
      aux2 =  (root_litter * froot_n2c) - aux2
      if(aawood .gt. 0.0) then
         aux3 = (cwd * awood_n2c) - aux3
      else
         aux3 = 0.0
      endif

      ! Check if it's nan
      if(isnan(aux1)) aux1 = 0.0
      if(isnan(aux2)) aux2 = 0.0
      if(isnan(aux3)) aux3 = 0.0
      ! ! Check if it's inf
      if(aux1 .eq. aux1 - 1) aux1 = 0.0
      if(aux2 .eq. aux2 - 1) aux2 = 0.0
      if(aux3 .eq. aux3 - 1) aux3 = 0.0

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

      if(aawood .gt. 0.0 .and. cwd .gt. 0.0) then
         litter_nutrient_ratios(3) = aux3 / cwd      ! N:C litter ratio g(N) g(C)-1
      else
         litter_nutrient_ratios(3) = 0.0
      endif


      ! P resobed
      aux1 = (leaf_litter * leaf_p2c) * rfrac_leaf
      aux2 = (root_litter * froot_p2c) * rfrac_froot
      if(aawood .gt. 0.0) then
         aux3 = (cwd * awood_p2c) * rfrac_wood
      else
         aux3 = 0.0
      endif
      ! ! Check if it's nan
      if(isnan(aux1)) aux1 = 0.0
      if(isnan(aux2)) aux2 = 0.0
      if(isnan(aux3)) aux3 = 0.0
      ! ! Check if it's inf
      if(aux1 .eq. aux1 - 1) aux1 = 0.0
      if(aux2 .eq. aux2 - 1) aux2 = 0.0
      if(aux3 .eq. aux3 - 1) aux3 = 0.0

      storage_out(3) = aux1 + aux2 + aux3

      aux1 =  (leaf_litter * leaf_p2c) - aux1    !nitrogen in litter
      aux2 =  (root_litter * froot_p2c) - aux2
      if(aawood .gt. 0.0) then
         aux3 = (cwd * awood_p2c) - aux3
      else
         aux3 = 0.0
      endif

      ! ! Check if it's nan
      if(isnan(aux1)) aux1 = 0.0
      if(isnan(aux2)) aux2 = 0.0
      if(isnan(aux3)) aux3 = 0.0
      ! ! Check if it's inf
      if(aux1 .eq. aux1 - 1) aux1 = 0.0
      if(aux2 .eq. aux2 - 1) aux2 = 0.0
      if(aux3 .eq. aux3 - 1) aux3 = 0.0

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

      if(aawood .gt. 0.0 .and. cwd .gt. 0.0) then
         litter_nutrient_ratios(6) = aux3 / cwd         ! P:C litter ratio g(N) g(C)-1
      else
         litter_nutrient_ratios(6) = 0.0
      endif


      ! END CALCULATIONS

      !    if(scl2_128 .lt. 0.0) scl2 = 0.0
      !    if(scf2_128 .lt. 0.0) scf2 = 0.0
      !    if(sca2_128 .lt. 0.0) sca2 = 0.0

      ! if(leaf_litter .lt. 0.0) leaf_litter = 0.0
      ! if(root_litter .lt. 0.0) root_litter = 0.0
      ! if(cwd .lt. 0.0) cwd = 0.0

      ! g m-2 to kg m-2
      !    scl2 = scl2 * 1e-3
      !    scf2 = scf2 * 1e-3
      !    sca2 = sca2 * 1e-3

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

      !    if(scl2_128 .lt. 0.0) scl2 = 0.0
      !    if(scf2_128 .lt. 0.0) scf2 = 0.0
      !    if(sca2_128 .lt. 0.0) sca2 = 0.0

      ! if(leaf_litter .lt. 0.0) leaf_litter = 0.0
      ! if(root_litter .lt. 0.0) root_litter = 0.0
      ! if(cwd .lt. 0.0) cwd = 0.0

      ! g m-2 to kg m-2

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
            if (k.eq.1) then
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
  !  subroutine carbon2 (tsoil,f5c,evap,laia,cl,cs)
  !     use types
  !     use global_par
  !     use photo_par
  !     implicit none
  !     !     Variables
  !     !     =========
  !     !     Inputs
  !     !     ------
  !     real(kind=r_4),intent(in) :: tsoil                !Mean monthly soil temperature (oC)
  !     real(kind=r_8),intent(in) :: f5c                  !Stress response to soil moisture (dimensionless)
  !     real(kind=r_8),intent(in) :: evap                 !Actual evapotranspiration (mm/day)
  !     real(kind=r_8),intent(in) :: laia
  !     !real(kind=r_4),intent(in) :: d_litter
  !     !     Outputs
  !     !     -------
  !     real(kind=r_4),dimension(2),intent(out) :: cl                   !Litter carbon (kgC/m2)
  !     real(kind=r_4),dimension(3),intent(out) :: cs                   !Soil carbon (kgC/m2)

  !     !     Internal
  !     !     --------
  !     real(kind=r_8) :: lf                   !Litterfall (kgC/m2)
  !     real(kind=r_8) :: f6                   !Litter decayment function
  !     real(kind=r_8) :: f7,zeta                   !Soil carbon storage function
  !     !
  !     !     Initialize
  !     !     ----------
  !     !
  !     lf  = 0.0
  !     f6  = 0.0
  !     f7  = 0.0
  !     cl  = 0.0
  !     cs  = 0.0

  !     !     Litter decayment function                                             !Controlled by annual evapotranspiration
  !     !     -------------------------
  !     f6 = 1.16*10.**(-1.4553+0.0014175*(evap*365.0))

  !     !     Soil carbon storage function                                          !Controlled by temperature
  !     !     ----------------------------
  !     f7 = p32**(p10*(tsoil-p11))

  !     !     Litterfall (kgC/m2)
  !     !     ------------------
  !     lf = p33 * laia ! + d_litter)

  !     !     Litter carbon (kgC/m2)
  !     !     ----------------------
  !     cl(1) = real(((lf/f6) * 0.75), r_4)
  !     cl(2) = real(((lf/f6) * 0.25), r_4)
  !     !     Soil carbon(kgC/m2)
  !     !     -------------------

  !     zeta = real(((p34 * (sum(cl)))/(p35*f7)) * f5c, kind=r_8)

  !     cs(1) = real(zeta * 0.60, r_4)
  !     cs(2) = real(zeta * 0.40, r_4)
  !     cs(3) = 0.0

  !     cl = cl * 1e3
  !     cs = cs * 1e3

  !  end subroutine carbon2

  !  !====================================================================
  !  !====================================================================

  !  ! Carbon decay implemented in JeDi and JSBACH - Pavlick et al. 2012

  !  function scarbon_decaiment(q10,tsoil,c,residence_time) result(decay)
  !     use types
  !     use global_par

  !     real(r_4),intent(in) :: q10              ! constant ~1.4
  !     real(r_4),intent(in) :: tsoil            ! Soil temperature °C
  !     real(r_4),intent(in) :: c                ! Carbon content per area g(C)m-2
  !     real(r_4),intent(in) :: residence_time   ! Pool turnover rate
  !     real(r_4) :: decay ! ML⁻²

  !     if(c .le. 0.0) then
  !        decay = 0.0
  !        return
  !     endif

  !     decay = (q10**((tsoil-20.0)/10.0)) * (c/residence_time)

  !  end function scarbon_decaiment

   !====================================================================
   !====================================================================

  !  subroutine carbon3(tsoil,leaf_l,cwd,root_l,lnr,cl,cs,cl_out,cs_out,hr)

  !     use types, only: r_4
  !     use global_par
  !     use soil_dec, only: nmin  => mineral_n_glob,&
  !          &  plab  => labile_p_glob,&
  !          &  l_run => last_run,&
  !          &  caux  => carbon_aux
  !     !implicit none
  !     integer(i_4),parameter :: pl=2,ps=3
  !     integer(i_4) :: index

  !     !     Variables
  !     !     =========
  !     !     Inputs
  !     !     ------
  !     real(r_4),intent(in) :: tsoil                 ! soil temperature (oC)
  !     real(r_4),intent(in) :: leaf_l                ! g(C)m⁻²
  !     real(r_4),intent(in) :: cwd
  !     real(r_4),intent(in) :: root_l
  !     real(r_4),dimension(6),intent(in) :: lnr      !g(Nutrient) g(C)⁻¹
  !     real(r_4),dimension(pl),intent(in) :: cl      !Litter carbon (gC/m2)
  !     real(r_4),dimension(ps),intent(in) :: cs      !Soil carbon (gC/m2)

  !     !     Outputs
  !     !     -------
  !     real(r_4),dimension(pl),intent(out) :: cl_out ! g(C)m⁻²
  !     real(r_4),dimension(ps),intent(out) :: cs_out
  !     real(r_4),intent(out) :: hr                   !Heterotrophic (microbial) respiration (gC/m2/day)

  !     !     Internal
  !     real(r_4),dimension(pl+ps) :: tr_c
  !     real(r_4),dimension(pl) :: pl_nitrogen = 0.0
  !     real(r_4),dimension(pl) :: pl_phosphorus = 0.0
  !     real(r_4),dimension(ps) :: ps_nitrogen = 0.0
  !     real(r_4),dimension(ps) :: ps_phosphorus = 0.0

  !     real(r_4) :: q10 = 1.4

  !     real(r_4) :: leaf_n2c
  !     real(r_4) :: froot_n2c
  !     real(r_4) :: wood_n2c
  !     real(r_4) :: leaf_p2c
  !     real(r_4) :: froot_p2c
  !     real(r_4) :: wood_p2c

  !     real(r_4) :: frac1,frac2
  !     real(r_4),dimension(5) :: het_resp, cdec
  !     real(r_4),dimension(5) :: aux_ratio_n, aux_ratio_p
  !     real(r_4),dimension(5) :: nutri_min_n, nutri_min_p
  !     real(r_4),parameter :: clit_atm = 0.77
  !     real(r_4),parameter :: cwd_atm = 0.2

  !     frac1 = 0.7
  !     frac2 = 1.0 - frac1

  !     ! Turnover Rates  == residence_time⁻¹ (years⁻¹)
  !     ! Coding residence time (years)
  !     tr_c(1) = 5.0    ! litter I   (1)
  !     tr_c(2) = 80.0   ! litter II  (2)
  !     tr_c(3) = 600.0   ! soil   I   (3)
  !     tr_c(4) = 2000.0  ! soil   II  (4)
  !     tr_c(5) = 1000.0 ! soil   III (5)

  !     ! find nutrient mass/area) : litter fluxes[ML⁻²] * litter nutrient ratios
  !     ! (lnr) [MM⁻¹]
  !     leaf_n2c = leaf_l * lnr(1) ! g(nutrient) m-2
  !     froot_n2c = root_l * lnr(2)
  !     wood_n2c = cwd * lnr(3)
  !     leaf_p2c = leaf_l * lnr(4)
  !     froot_p2c = root_l * lnr(5)
  !     wood_p2c = cwd * lnr(6)

  !     ! FIRST OF ALL calculate dacay from pools
  !     ! CARBON DECAY
  !     do index = 1,4
  !        if(index .lt. 3) then
  !           ! FOR THE 2 LITTER POOLS
  !           cdec(index) = scarbon_decaiment(q10,tsoil,cl(index),tr_c(index))
  !        else
  !           ! FOR THE 3 CARBON POOLS
  !           cdec(index) = scarbon_decaiment(q10,tsoil,cs(index-2),tr_c(index))
  !        endif
  !     enddo
  !     cdec(5) = 0.0

  !     ! partitioning material coming from vegetation
  !     ! filling C (litter I-II; soil I) pools with incoming material
  !     !LITTER I & II
  !     cl_out(1) = (cl(1) - cdec(1)) + (frac1 * leaf_l) + (frac1 * root_l)
  !     cl_out(2) = (cl(2) - cdec(2)) + (frac2 * leaf_l) + (frac2 * root_l)
  !     cl_out(2) = cl_out(2) + (cwd * frac2)
  !     !SOIL Ie
  !     cs_out(1) = (cs(1) - cdec(3)) + (frac1 * cwd)
  !     cs_out(2) = ((1 - clit_atm) * (cdec(1) + cdec(2))) + &
  !          & ((1 - cwd_atm) * cdec(3)) - cdec(4)
  !     cs_out(3) = 0.0 ! cs_out(3) is empty
  !     ! print *, 'cl after partitioning->', cl
  !     ! print *, 'cs after partitioning->', cs

  !     ! Calculate the nutrint contents in each litterI-II/soilI pool
  !     pl_nitrogen(1)   = (leaf_n2c * frac1)+(froot_n2c * frac1)  ! g(N)m-2
  !     pl_nitrogen(2)   = (leaf_n2c * frac2)+(froot_n2c * frac2)+(wood_n2c *frac2)
  !     ps_nitrogen(1)   = wood_n2c * frac1
  !     pl_phosphorus(1) = (leaf_p2c * frac1)+(froot_p2c * frac1)
  !     pl_phosphorus(2) = (leaf_p2c * frac2)+(froot_p2c * frac2)+(wood_p2c*frac2)
  !     ps_phosphorus(1) = wood_p2c * frac1
  !     ! System Memory !WARNING! carbon aux is a global variable
  !     ! I do not know why the hell i did it
  !     ps_nitrogen(2) = caux(1)
  !     ps_nitrogen(3) = 0.0
  !     ps_phosphorus(2) = caux(3)
  !     ps_phosphorus(3) = 0.0
  !     ! print *, 'n in litter', pl_nitrogen
  !     ! print *, 'p in litter', pl_phosphorus
  !     ! print *, 'n in soil', ps_nitrogen
  !     ! print *, 'p in soil', ps_phosphorus

  !     !HRESP
  !     het_resp(1) = clit_atm * cdec(1)
  !     het_resp(2) = clit_atm * cdec(2)
  !     het_resp(3) = cwd_atm * cdec(3)
  !     het_resp(4) = cdec(4)
  !     het_resp(5) = 0.0

  !     ! Mineralized nutrients
  !     do index=1,2
  !        if(cl_out(index) .gt. 0.0) then
  !           aux_ratio_n(index) = pl_nitrogen(index)/cl_out(index)   ! g(N)g(C)-1
  !           aux_ratio_p(index) = pl_phosphorus(index)/cl_out(index) ! g(P)g(C)-1
  !        else
  !           aux_ratio_n(index) = 0.0
  !           aux_ratio_p(index) = 0.0
  !        endif
  !     enddo

  !     do index=3,4
  !        if(cs_out(index-2) .gt. 0.0) then
  !           aux_ratio_n(index) = ps_nitrogen(index-2)/cs_out(index-2)   ! g(N)g(C)-1
  !           aux_ratio_p(index) = ps_phosphorus(index-2)/cs_out(index-2) ! g(P)g(C)-1
  !        else
  !           aux_ratio_n(index) = 0.0
  !           aux_ratio_p(index) = 0.0
  !        endif
  !     enddo
  !     aux_ratio_p(5) = 0.0

  !     ! print *, 'ratiop',aux_ratio_p
  !     ! print *, 'ration',aux_ratio_n

  !     do index=1,4
  !        nutri_min_n(index) = het_resp(index)*aux_ratio_n(index)
  !        nutri_min_p(index) = het_resp(index)*aux_ratio_p(index)
  !     enddo
  !     nutri_min_n(5) = 0.0
  !     nutri_min_p(5) = 0.0

  !     do index=1,2
  !        pl_nitrogen(index) = pl_nitrogen(index) - nutri_min_n(index)
  !        pl_phosphorus(index) = pl_phosphorus(index) - nutri_min_p(index)
  !     enddo

  !     do index=1,3
  !        ps_nitrogen(index) = ps_nitrogen(index) - nutri_min_n(index + 2)
  !        ps_phosphorus(index) = ps_phosphorus(index) - nutri_min_p(index + 2)
  !     enddo


  !     ! calculate final state

  !     ! SOIL II and III nutrient pools
  !     caux(1) =  ps_nitrogen(2)
  !     caux(2) =  0.0
  !     caux(3) =  ps_phosphorus(2)
  !     caux(4) =  0.0

  !     if(l_run) then
  !        continue
  !        ! fill soil_dec_variables
  !     endif

  !     ! FEATURES TO BE IMPLEMENTED

  !     ! P weathering +
  !     ! P biochemical mineralization +
  !     ! P release "de-sorption" +
  !     ! P immobilization
  !     ! P leaching -
  !     ! P sorption I & II -
  !     ! P occlusion -

  !     ! N deposition +
  !     ! N fixation +
  !     ! N immobilization -
  !     ! N leaching -
  !     ! N degassing-volatilization -

  !     nmin = nmin + (sum(nutri_min_n,&
  !          &mask=.not. isnan(nutri_min_n)) * 1e-3) ! - kg m⁻²
  !     plab = plab + (sum(nutri_min_p,&
  !          &mask=.not. isnan(nutri_min_p)) * 1e-3) ! - kg m⁻²

  !     hr = sum(het_resp)

  !     ! print *, 'c AUX', carbon_aux
  !     ! print *, 'n mineral', nutri_min_n
  !     ! print *, 'p labila', nutri_min_p

  !     ! print *, 'carbon3 _pools'
  !     ! print *, 'plab ---.>',plab
  !     ! print *, 'nmin ---.>',nmin
  !     ! print *, 'hr ---.>',hr


  !  end subroutine carbon3

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



   !   subroutine pft_area_frac(cleaf1, cfroot1, cawood1, ocp_coeffs, ocp_wood)
   !     use types, only: l_1, i_4, r_4, r_8
   !     use global_par, only: npls, sapwood,cmin
   !     !implicit none

   !     integer(i_4),parameter :: npft = npls ! plss futuramente serao

   !     real(r_8),dimension(npft),intent( in) :: cleaf1, cfroot1, cawood1
   !     real(r_8),dimension(npft),intent(out) :: ocp_coeffs
   !     logical(l_1),dimension(npft),intent(out) :: ocp_wood

   !     real(r_8),dimension(npft) :: cleaf, cawood, cfroot
   !     real(r_8),dimension(npft) :: total_biomass_pft=0.0
   !     real(r_8),dimension(npft) :: total_w_pft = 0.0
   !     integer(i_4) :: p,i
   !     integer(i_4),dimension(1) :: max_index
   !     real(r_8) :: total_biomass = 0.0
   !     real(r_8) :: total_wood = 0.0
   !     integer(i_4) :: five_percent = 0

   !     ! ! Scalars
   !     ! total_biomass = 0.0
   !     ! total_wood = 0.0

   !     ! # Arrays
   !     cleaf = cleaf1
   !     cfroot = cfroot1
   !     cawood = cawood1
   !     ! total_w_pft = 0.0
   !     ! total_biomass_pft = 0.0
   !     ocp_coeffs = 0.0
   !     ocp_wood = .false.

   !     do p = 1,npft
   !        ocp_coeffs(p) = 0.0
   !        ocp_wood(p) = .false.
   !     enddo

   !     ! check for nan in cleaf cawood cfroot
   !     do p = 1,npft
   !        if(isnan(cleaf(p))) cleaf(p) = 0.0
   !        if(isnan(cfroot(p))) cfroot(p) = 0.0
   !        if(isnan(cawood(p))) cawood(p) = 0.0

   !        if(cleaf(p) .eq. cleaf(p) - 1) cleaf(p) = 0.0
   !        if(cfroot(p) .eq. cfroot(p) - 1) cfroot(p) = 0.0
   !        if(cawood(p) .eq. cawood(p) - 1) cawood(p) = 0.0
   !     enddo

   !     do p = 1,npft
   !       ! ! ! check for PLS availabilityo
   !       ! !  if(((cleaf(p) .lt. cmin) .and. (cfroot(p) .lt. cmin))) then
   !       !     total_biomass_pft(p) = 0.0
   !       !     total_w_pft(p) = 0.0
   !       !  !    goto 795
   !       !  ! endif
   !        ! Then Calculates
   !        total_biomass_pft(p) = cleaf(p) + cfroot(p) + cawood(p) ! only sapwood
   !        total_w_pft(p) = cawood(p)

   ! ! 795    continue
   ! !        total_biomass = total_biomass + 0.0
   ! !        total_wood = total_wood 0.0
   ! !     enddo



   !     !     grid cell occupation coefficients
   !     if(total_biomass .gt. 0.0) then
   !        do p = 1,npft
   !           ocp_coeffs(p) = total_biomass_pft(p) / total_biomass
   !           if(ocp_coeffs(p) .lt. 0.0) ocp_coeffs(p) = 0.0
   !           if(isnan(ocp_coeffs(p))) ocp_coeffs(p) = 0.0
   !        enddo
   !     else
   !        do p = 1,npft
   !           ocp_coeffs(p) = 0.0
   !        enddo
   !     endif

   !     !     gridcell pft ligth limitation by wood content
   !     five_percent = nint((real(npft) * 0.05),kind=i_4)
   !     if(five_percent .eq. 0) five_percent = 1
   !     if(five_percent .eq. 1) then
   !        if(total_wood .gt. 0.0) then
   !           max_index = maxloc(total_w_pft)
   !           i = max_index(1)
   !           ocp_wood(i) = .true.
   !        endif
   !     else
   !        do p = 1,five_percent
   !           if(total_wood .gt. 0.0) then
   !              max_index = maxloc(total_w_pft)
   !              i = max_index(1)
   !              total_w_pft(i) = 0.0
   !              ocp_wood(i) = .true.
   !           endif
   !        enddo
   !     endif

   !   return
   !   end subroutine pft_area_frac


   !====================================================================
   ! ALGUMAS FUNÇÕES QUE SÃO UTEIS AS VEZES (QUASE NUNCA)
   !====================================================================

   subroutine pft_par(dt)
      ! DEPRECATED
      use types, only: r_4
      use global_par, only: ntraits, npls
      !implicit none


      real(r_4), dimension(npls,ntraits),intent(out) :: dt

      !     dt1 = g1       kPa ** (1/2)
      !     dt2 = vcmax    mol m-2 s-1
      !     dt3 = tleaf    years
      !     dt4 = twood    years
      !     dt5 = tfroot   years
      !     dt6 = aleaf    unitless
      !     dt7 = awood    unitless
      !     dt8 = aroot    unitless
      !     dt9 = is C4 [1 = yes, 0 = no]
      !     dt10 = leaf N:C
      !     dt11 = wood N:C
      !     dt12 = froot N:C
      !     dt13 = leaf P:C
      !     dt14 = awood P:C
      !     dt15 = froot P:C

      open(45,file='pls.bin',status='old',&
           &form='unformatted',access='direct',recl=4*npls*ntraits)
      read(45,rec=1) dt
      close(45)
      return
   end subroutine pft_par

   !=================================================================
   !=================================================================

   subroutine ascii2bin(file_in, file_out, nx1, ny1)
      ! DEPRECATED
      use types
      !implicit none

      character*100, intent(in) :: file_in, file_out
      integer(i_4),intent(in) :: nx1, ny1

      integer(i_4) :: i, j

      real(r_4),allocatable,dimension(:,:) :: arr_in

      allocate(arr_in(nx1,ny1))

      open (unit=11,file=file_in,status='old',form='formatted',access='sequential',&
           action='read')


      open (unit=21,file=file_out,status='unknown',&
           form='unformatted',access='direct',recl=nx1*ny1*r_4)


      do j = 1, ny1 ! for each line do
         read(11,*) (arr_in(i,j), i=1,nx1) ! read all elements in line j (implicitly looping)
         !write(*,*) arr_in(:,j)
      end do

      write(21,rec=1) arr_in
      close(11)
      close(21)

   end subroutine ascii2bin

   !=================================================================
   !=================================================================

   function leap(year) result(is_leap)
      ! Why this is here?
      use types
      !implicit none

      integer(i_4),intent(in) :: year
      logical(l_1) :: is_leap

      logical(l_1) :: by4, by100, by400

      by4 = (mod(year,4) .eq. 0)
      by100 = (mod(year,100) .eq. 0)
      by400 = (mod(year,400) .eq. 0)

      is_leap = by4 .and. (by400 .or. (.not. by100))

   end function leap

   !=================================================================
   !=================================================================

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
