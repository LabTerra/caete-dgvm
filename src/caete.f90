! Copyright 2017- LabTerra

!     This program is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Fundation, either version 3 of the License, or
!     (at your option) any later version.)

!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.

!     You should have received a copy of the GNU General Public License
!     along with this program.  If not, see <http://www.gnu.org/licenses/>.

module caete
   implicit none
   private

   public ::  caete_dyn


   contains

   subroutine caete_dyn(x, y, run, dt, w0, g0, s0, dcl, dca, dcf, prec,&
        & temp, p0, par, rhs, cleaf_ini, cawood_ini, cfroot_ini, emaxm, tsoil,&
        & photo_cwm, aresp_cwm, npp_cwm, lai_cwm, hr_cwm, rcm_cwm, f51,&
        & runom_cwm, evapm_cwm, wsoil_cwm, rm_cwm, rg_cwm, cleaf_cwm,&
        & cawood_cwm, cfroot_cwm, grid_area, wue, cue, cdef, wfim, gfim, sfim,&
        & dl_final, dw_final, dr_final, clf, caf, cff, vcmax, specific_la,&
        & nupt, pupt,litter_l, cwd, litter_fr, lnr, storage_pool, soil_carbon,&
        & snr, avail_p, inorg_n, inorg_p, sorbed_p)

      use types
      use utils, only: gpid   => process_id
      use global_par
      use photo, only: cwm4, cwm8, cwm_soil
      use water, only : soil_temp, soil_temp_sub
      use budget, only : daily_budget
      use soil_dec
      use soil_pools

      !     --------------------------I N P U T S----------------------------

      !boundary conditions - initial state
      integer(i_4), intent(in) :: x, y, run  ! gridcell ID
      ! integer(i_4), intent(in) :: ndays      ! Number of simulated days
      real(r_4),dimension(ntraits, npls),intent(in) :: dt ! trait table
      real(r_4),dimension(npls),intent(in) :: w0 ! Initial soil water kg m-2
      real(r_4),dimension(npls),intent(in) :: g0 ! Initial soil ice
      real(r_4),dimension(npls),intent(in) :: s0 ! Initial soil snow

      ! , csoil_init, snr_init,&
      !  & in_p_init, in_n_init, av_p_init, so_p_init
      ! real(r_4),dimension(4),intent(in) :: csoil_init
      ! real(r_4),dimension(8),intent(in) :: snr_init
      ! real(r_4),intent(in) :: in_p_init
      ! real(r_4),intent(in) :: in_n_init
      ! real(r_4),intent(in) :: av_p_init
      ! real(r_4),intent(in) :: so_p_init

      ! initial value for betaleaf (daily delta cmass) Gained C - Used to calculate Growth Respiration
      real(r_8),dimension(npls),intent(in) :: dcl ! Leaf - deltaC(day, day-1) kg m-2
      real(r_8),dimension(npls),intent(in) :: dca ! wood tissues- deltaC(day, day-1) kg m-2
      real(r_8),dimension(npls),intent(in) :: dcf ! fine root - deltaC(day, day-1) kg m-2

      ! Meteorological variables - INPUT DATA
      real(r_4),dimension(nt1),intent(in) :: p0         ! entry: Pa; convert to hPa - Atmospheric pressure (mb)
      real(r_4),dimension(nt1),intent(in) :: prec       ! entry: Kg m-2 s-1; convert to mm month-1 - Precipitation (mm/month)
      real(r_4),dimension(nt1),intent(in) :: temp       ! entry: K; convert - Temperature (oC)
      real(r_4),dimension(nt1),intent(in) :: par        ! entry: RSDS shortwave downward rad. -W m-2; convert - IPAR (Ein/m2/s) molþ m-2 s-1
      real(r_4),dimension(nt1),intent(in) :: rhs        ! entry: %; tranform to ratio (0-1) Relative humidity

      real(r_8),dimension(npls),intent(in) :: cleaf_ini    ! Initial carbon content in leaves (kg m-2)
      real(r_8),dimension(npls),intent(in) :: cawood_ini   ! Initial carbon content in aboveground wood (kg m-2)
      real(r_8),dimension(npls),intent(in) :: cfroot_ini   ! Initial carbon content in fineroots (kg m-2)

      !     -----------------------------E N D-------------------------------

      !     -------------------------O U T P U T S---------------------------

      real(r_4),dimension(nt1),   intent(out) :: emaxm        ! 1   ! Max.evapotranspiration (kg m-2 day-1)
      real(r_4),dimension(nt1),   intent(out) :: tsoil        ! 2   ! soil temperature °C
      real(r_4),dimension(nt1),   intent(out) :: photo_cwm    ! 3   ! daily photosynthesis   (kgC m-2 year-1)
      real(r_4),dimension(nt1),   intent(out) :: aresp_cwm    ! 4   ! daily autotrophic res  (kgC m-2 year-1)
      real(r_4),dimension(nt1),   intent(out) :: npp_cwm      ! 5   ! daily net primary produ (kgC m-2 year-1)
      real(r_4),dimension(nt1),   intent(out) :: lai_cwm      ! 6   ! daily leaf area index m2 m-2
      real(r_4),dimension(nt1),   intent(out) :: hr_cwm       ! 7   ! daily het resp  (kgC/m2)
      real(r_8),dimension(nt1),   intent(out) :: rcm_cwm      ! 8  ! leaf resistence s m-1
      real(r_8),dimension(nt1),   intent(out) :: f51          ! 9  ! Water stress modifier (dimensionless)
      real(r_8),dimension(nt1),   intent(out) :: runom_cwm    ! 10  ! Runoff kg m-2 day-1
      real(r_8),dimension(nt1),   intent(out) :: evapm_cwm    ! 11  ! Actual evapotranspiration kg m-2 day-1
      real(r_8),dimension(nt1),   intent(out) :: wsoil_cwm    ! 12  ! Soil moisture (kg m-2) check in transition (run_dyn)
      real(r_8),dimension(nt1),   intent(out) :: rm_cwm       ! 13  ! Plant (autotrophic) Maintenance respiration
      real(r_8),dimension(nt1),   intent(out) :: rg_cwm       ! 14  ! Plant (autotrophic) Growth respiration
      real(r_8),dimension(nt1),   intent(out) :: cleaf_cwm    ! 15  ! leaf biomass (KgC/m2)
      real(r_8),dimension(nt1),   intent(out) :: cawood_cwm   ! 16  ! aboveground wood biomass (KgC/m2)
      real(r_8),dimension(nt1),   intent(out) :: cfroot_cwm   ! 17  ! fine root biomass (KgC/m2)
      real(r_8),dimension(npls),  intent(out) :: grid_area    ! 18  ! gridcell area fraction of pfts! ratio (0-1)
      real(r_8),dimension(nt1),   intent(out) :: wue          ! 19  ! 1 - molCO2 m-2 s-1 (molH2O m-2 s-1)-1;
      real(r_8),dimension(nt1),   intent(out) :: cue          ! 20  ! 2 - npp/gpp;
      real(r_8),dimension(nt1),   intent(out) :: cdef         ! 21  ! 3 - g m-2 day-1
      real(r_4),dimension(npls),  intent(out) :: wfim         ! 22  ! Kg m-2
      real(r_4),dimension(npls),  intent(out) :: gfim         ! 23  ! Kg m-2
      real(r_4),dimension(npls),  intent(out) :: sfim         ! 24  ! Kg m-2  final day water pools (snow; water; ice)

      real(r_8),dimension(npls),  intent(out) :: dl_final     ! 25  ! delta(now, last day pool size)
      real(r_8),dimension(npls),  intent(out) :: dw_final     ! 26  !
      real(r_8),dimension(npls),  intent(out) :: dr_final     ! 27  ! KgC m-2
      real(r_8),dimension(npls),  intent(out) :: clf          ! 28  ! final carbon pools (cveg) for each pft (not scaled to cell area)
      real(r_8),dimension(npls),  intent(out) :: caf          ! 29  ! KgC m-2
      real(r_8),dimension(npls),  intent(out) :: cff          ! 30

      real(r_4),dimension(nt1),   intent(out) :: vcmax        ! 31  ! mol m⁻² s⁻¹
      real(r_4),dimension(nt1),   intent(out) :: specific_la  ! 32  ! m²g⁻¹
      real(r_4),dimension(nt1),   intent(out) :: nupt         ! 33  ! n plant uptake - CWM  TODO units
      real(r_4),dimension(nt1),   intent(out) :: pupt         ! 34  ! p plant uptake - CWM
      real(r_4),dimension(nt1),   intent(out) :: litter_l     ! 35  ! CWM fluxes from vegetation to soil g m⁻² day⁻¹
      real(r_4),dimension(nt1),   intent(out) :: cwd          ! 36  ! CWM
      real(r_4),dimension(nt1),   intent(out) :: litter_fr    ! 37  ! CWM

      ! Litter Nutrient Ratisos :: variables(6)         [(lln2c),(rln2c),(cwdn2c),(llp2c),(rlp2c),(cwdp2c)]
      real(r_4),dimension(6,nt1), intent(out) :: lnr          ! 38 ! Litter nutrient ratios

      real(r_4),dimension(3,nt1), intent(out) :: storage_pool ! 39 ! Rapid Auxiliary daily storage pool for carbon and nutrients g m⁻²
      real(r_4),dimension(4,nt1), intent(out) :: soil_carbon  ! 40 ! CSOIL_POOL carbon gC m-2

     ! Soil Nutrient Ratios :: variables(8)            [(l1n2c),(l2n2c),(c1dn2c),(c2n2c),(l1p2c),(l2p2c),(c1p2c),(c2p2c)]
      real(r_4),dimension(8,nt1), intent(out) :: snr          ! 41 ! Soil nutrient ratios g/g
      real(r_4),dimension(nt1),   intent(out) :: avail_p      ! 42
      real(r_4),dimension(nt1),   intent(out) :: inorg_n      ! 43
      real(r_4),dimension(nt1),   intent(out) :: inorg_p      ! 44
      real(r_4),dimension(nt1),   intent(out) :: sorbed_p     ! 45

      !  c     --------------------------------E N D----------------------------

      !  c     ------------------------- internal variables---------------------
      integer(i_4) :: k, index
      real(r_8),dimension(nt1) :: gsoil    !Soil ice kg m-2
      real(r_8),dimension(nt1) :: ssoil    !Soil snow kg m-2
      real(r_8),dimension(nt1) :: snowm    !Snowmelt kg m-2

      real(r_8),dimension(npls) :: vcmax_com
      real(r_8),dimension(npls) :: specific_la_com
      real(r_8),dimension(npls) :: nupt_com
      real(r_8),dimension(npls) :: pupt_com
      real(r_8),dimension(npls) :: litter_l_com
      real(r_8),dimension(npls) :: cwd_com
      real(r_8),dimension(npls) :: litter_fr_com
      real(r_8),dimension(6,npls) :: lnr_com

      real(r_4),dimension(npls) :: sini  !kg m-2 INPUT TO daily_budget
      real(r_4),dimension(npls) :: wini  !kg m-2 ! Change in daily loop
      real(r_4),dimension(npls) :: gini  !kg m-2

      ! FEED AN INOUT VARIABLE IN daily_budget
      real(r_8),dimension(3,npls) :: storage_pool_com ! FEED AN INOUT VARIABLE IN daily_budget
      real(r_8),dimension(npls) :: dl  ! delta(now, last day pool size)
      real(r_8),dimension(npls) :: dr  ! kg m-2 INPUT TO daily_budget
      real(r_8),dimension(npls) :: dw  ! Change in daily loop
      real(r_8),dimension(npls) :: cleaf1_pft   !kg m-2 ! FEED AN INOUT VARIABLE to daily_budget
      real(r_8),dimension(npls) :: cawood1_pft  !kg m-2 !
      real(r_8),dimension(npls) :: cfroot1_pft  !kg m-2

      ! Recv the outputs from daily_budget
      real(r_4)                 :: ep_daily                !potential evapotranspiration(PET)(mm day-1)
      real(r_8),dimension(npls) :: runoff_com            ! mm day-1
      real(r_8),dimension(npls) :: gpp_com               ! Kg m-2 year-1
      real(r_8),dimension(npls) :: snowmelt_com         ! mm day-1
      real(r_8),dimension(npls) :: canopy_res_com ! s m-1
      real(r_8),dimension(npls) :: f5_com
      real(r_8),dimension(npls) :: npp_com
      real(r_8),dimension(npls) :: lai_com
      real(r_8),dimension(npls) :: ar_com
      real(r_8),dimension(npls) :: e_com
      real(r_8),dimension(npls) :: rm_com
      real(r_8),dimension(npls) :: rg_com
      real(r_8),dimension(npls) :: cue_com
      real(r_8),dimension(npls) :: c_def_com
      real(r_8),dimension(npls) :: cleafcom
      real(r_8),dimension(npls) :: cawoodcom
      real(r_8),dimension(npls) :: cfrootcom
      real(r_8),dimension(npls) :: gridocpcom
      real(r_4),dimension(npls) :: hr_com
      real(r_8),dimension(npls) :: wue_com

      real(r_8),dimension(4,npls) :: soilcarbon_com
      real(r_4),dimension(8,npls) :: snr_com
      real(r_8),dimension(npls)   :: in_p_com
      real(r_8),dimension(npls)   :: in_n_com
      real(r_8),dimension(npls)   :: av_p_com
      real(r_8),dimension(npls)   :: so_p_com

      ! Auxiliary to loop transition
      real(r_4),dimension(4) :: csoil_in
      real(r_4),dimension(8) :: snr_in
      real(r_4) :: in_p_in, in_n_in, av_p_in, so_p_in ! Auxiliary variables
      real(r_8),dimension(npls) :: grd

      real(r_4) :: pr,spre,ta,td,ipar,ru,negative_n, negative_p
      integer(i_4) :: p
      real(r_4) :: t1ww
      real(r_4) :: t2ww
      real(r_4) :: t3ww

      ! Next are auxiliary to tests
      integer(i_4),dimension(npls) :: gridocpcom_int
      integer(i_4) :: mypid_int, ls, npq, spnp
      character(len=29) :: filename
      character(len=7) :: mypid
      character(len=3) :: str_x, str_y
      character(len=9), parameter :: fname = 'pools_din'
      character(len=4), parameter :: fext = '.txt'
      !  END OF DECLARATION

      ! ### MODEL INITIALIZATION ---------------------

      ! Create a textfile for debug
      if(text_ts) then
         ! open log file
         mypid_int = gpid()
         write (mypid, '(I7)') mypid_int
         write (str_x, '(I3)') x
         write (str_y, '(I3)') y
         filename =  adjustl(fname // trim(adjustl(mypid)) //'-' // trim(adjustl(str_x))&
            & // '-' // trim(adjustl(str_y)) // fext)

         open (unit=12345,file=trim(filename),status='unknown',&
              &form='formatted',position='append', action='write')
      endif

      ! Carbon pools initial values
      cleaf1_pft  =  cleaf_ini    ! Carbon pools
      cawood1_pft = cawood_ini
      cfroot1_pft = cfroot_ini

      ! 'Growth' pools (to be used by growth respiration)
      dl = dcl                   ! Daily Delta in carbon pools
      dr = dcf
      dw = dca

      ! Soil Waters pools
      wini  = w0  !Soil moisture_initial condition (mm)
      gini  = g0  !Soil ice_initial condition (mm)
      sini  = s0  !Overland snow_initial condition (mm)

      ! Soil organic matter pools
      in_p_in  = sp_in_p
      in_n_in  = sp_available_n
      av_p_in  = sp_available_p
      so_p_in  = sp_so_p
      csoil_in = sp_csoil
      snr_in   = sp_snr

      !82 columns-----------------------------------------------------------------

      ! Initialize variables OUTPUTS
      emaxm              = 0.0     ! Maximum evapotranspiration
      tsoil              = 0.0
      photo_cwm          = 0.0     ! daily photosynthesis (kgC/m2/y)
      aresp_cwm          = 0.0     ! daily autotrophic respiration (kgC/m2/y)
      npp_cwm            = 0.0     ! daily net primary productivity (sum of PLSs) (kgC/m2/y)
      lai_cwm            = 0.0     ! daily leaf area index m2 m-2n
      hr_cwm             = 0.0     ! daily heterotrophic respiration (kgC/m2/y)
      rcm_cwm            = 0.0     ! molH2O m-2 s-1
      f51                = 0.0     ! dimensionless
      gsoil              = 0.0     ! Soil ice mm
      ssoil              = 0.0     ! Soil snow mm
      runom_cwm          = 0.0     ! Runoff mm
      evapm_cwm          = 0.0     ! Actual evapotranspiration mm/day
      wsoil_cwm          = 0.0     ! Soil moisture (mm)
      rm_cwm             = 0.0     ! Maintenance respiration
      rg_cwm             = 0.0
      cleaf_cwm          = 0.0     ! leaf biomass (KgC/m2)
      cawood_cwm         = 0.0     ! aboveground wood biomass (KgC/m2)
      cfroot_cwm         = 0.0     ! fine root biomass (KgC/m2)
      grid_area          = 0.0     ! occupation coefficient
      wue                = 0.0     ! Water use efficiency (Medlyn et al. 2011 Reconciling...)
      cue                = 0.0     ! Carbon Use efficiency (NPP/GPP)
      cdef               = 0.0

      !  ### End model initialization -------------------------------------------------

         ! Daily budget
         ! ====================
         ! Variables that I will fill with data in call daily_budget
      soilcarbon_com          = 0.0     ! Carbon in soil g C m-2 (l1, l2, s1, s2)
      snr_com                 = 0.0
      in_p_com                = 0.0
      in_n_com                = 0.0
      av_p_com                = 0.0
      so_p_com                = 0.0
      cleafcom                = 0.0
      cawoodcom               = 0.0
      cfrootcom               = 0.0
      gridocpcom              = 0.0
      hr_com                  = 0.0


      !     ======================
      !     START TIME INTEGRATION
      !     ======================
      do k = 1,nt1
         spnp = 0
         if (k .gt. 1) then
            tsoil(k) = soil_temp(tsoil(k-1), temp(k) - 273.15)
         else
            call soil_temp_sub(temp(1:1095) - 273.15,tsoil(k))
         endif

                  ! Daily budget
         ! ====================
         ! Variables that I will fill with data in call daily_budget
         ep_daily                = 0.0

         ! dimension(NPLS) - OUTPUTS
         !wfim                    = 0.0
         !gfim                    = 0.0
         !sfim                    = 0.0
         !dl_final                = 0.0
         !dr_final                = 0.0
         !dw_final                = 0.0
         !clf                     = 0.0
         !caf                     = 0.0
         !cff                     = 0.0
         ! INTERNAL - CATCH outputs from COMMUNITY after budget
         ! This section clean the values for each day
         ! dimension(NPLS)
         snowmelt_com            = 0.0
         runoff_com              = 0.0
         e_com                   = 0.0
         gpp_com                 = 0.0
         ar_com                  = 0.0
         npp_com                 = 0.0
         lai_com                 = 0.0
         f5_com                  = 0.0
         canopy_res_com          = 0.0
         rm_com                  = 0.0
         rg_com                  = 0.0
         wue_com                 = 0.0
         cue_com                 = 0.0
         c_def_com               = 0.0

         soilcarbon_com          = 0.0     ! Carbon in soil g C m-2 (l1, l2, s1, s2)
         snr_com                 = 0.0
         in_p_com                = 0.0
         in_n_com                = 0.0
         av_p_com                = 0.0
         so_p_com                = 0.0
         cleafcom                = 0.0
         cawoodcom               = 0.0
         cfrootcom               = 0.0
         gridocpcom              = 0.0
         hr_com                  = 0.0

         ! Setting input variables and Converting units
         td   = tsoil(k)
         spre = p0(k) * 0.01            ! transforamando de Pascal pra mbar (hPa)
         ta   = temp(k) - 273.15        ! K to °C
         pr   = prec(k) * 86400.0         !2.62974e+06 ! kg m-2 s-1 to mm month-1
         ipar = (0.5 * par(k)) / 2.18e5 ! W m-2 to mol m-2 s-1 ! 0.5 converts RSDS to PAR
         ru   = rhs(k) / 100.0          ! Relative humidity


         ! if(k .lt. 30000) then
            ! if (in_n_in .gt. sp_available_n) in_n_in = sp_available_n
            ! if (av_p_in .gt. sp_available_p) av_p_in = sp_available_p
         ! endif

         call daily_budget(dt, wini, gini, sini, td, ta, pr, spre, ipar, ru &
              &, in_p_in, in_n_in, av_p_in, so_p_in, storage_pool_com, cleaf1_pft &
              &, cawood1_pft, cfroot1_pft, dl, dw, dr, csoil_in, snr_in, wfim, gfim, sfim &
              &, snowmelt_com, runoff_com, e_com, ep_daily, gpp_com, ar_com, npp_com &
              &, lai_com, canopy_res_com, f5_com, rm_com, rg_com, cleafcom, cawoodcom &
              &, cfrootcom, gridocpcom, wue_com, cue_com, c_def_com, vcmax_com &
              &, specific_la_com, soilcarbon_com, in_p_com, in_n_com, av_p_com, so_p_com &
              &, nupt_com, pupt_com, litter_l_com, cwd_com &
              &, litter_fr_com, hr_com, lnr_com, snr_com)


         !82 columns-------------------------------------------------------------
         !print *, cawoodcom, 'cw'
         !print *, gridocpcom, 'grd'
         grd            = gridocpcom
         emaxm(k)       = ep_daily
         gsoil(k)     = sum(gfim * grd,    mask= .not. isnan(gfim))
         ssoil(k)     = sum(sfim * grd,    mask= .not. isnan(sfim))
         snowm(k)     = sum(snowmelt_com * grd,    mask= .not. isnan(snowmelt_com))

         wsoil_cwm(k) = real(sum(wfim * grd,    mask= .not.  isnan(wfim)),r_4)
         runom_cwm(k) = real(sum(runoff_com * grd,    mask= .not.  isnan(runoff_com)),r_4)
         evapm_cwm(k) = real(sum(e_com * grd,    mask= .not.  isnan(e_com)),r_4)
         f51(k)       = real(sum(f5_com * grd,    mask= .not. isnan(f5_com)),r_4)

         rcm_cwm(k)   = real(sum(canopy_res_com * grd,   mask= .not. isnan(canopy_res_com)),r_4)
         lai_cwm(k)   = real(sum(lai_com* grd,   mask= .not. isnan(lai_com)),r_4)
         photo_cwm(k) = real(sum(gpp_com * grd,   mask= .not. isnan(gpp_com)),r_4)
         aresp_cwm(k) = real(sum(ar_com * grd,   mask= .not. isnan(ar_com)),r_4)
         npp_cwm(k)   = real(sum(npp_com* grd,   mask= .not. isnan(npp_com)),r_4)
         rm_cwm(k)    = real(sum(rm_com * grd,   mask= .not. isnan(rm_com)),r_4)
         rg_cwm(k)    = real(sum(rg_com * grd,   mask= .not. isnan(rg_com)),r_4)
         wue(k)       = real(sum(wue_com* grd,   mask= .not. isnan(wue_com)),r_4)
         cue(k)       = real(sum(cue_com* grd,   mask= .not. isnan(cue_com)),r_4)
         cdef(k)      = real(sum(c_def_com* grd, mask= .not. isnan(c_def_com)),r_4)

         vcmax(k)  = real(sum(vcmax_com * grd,  mask= .not. isnan(vcmax_com)),r_4)

         specific_la(k) = real(sum(specific_la_com  * grd,&
              & mask=.not. isnan(specific_la_com)),r_4)

         nupt(k) = real(sum(nupt_com * grd, mask= .not. isnan(nupt_com)),r_4)
         pupt(k) = real(sum(pupt_com * grd, mask= .not. isnan(pupt_com)),r_4)

         litter_l(k) = real(sum(litter_l_com * grd,&
              & mask=.not. isnan(litter_l_com)),r_4)
         cwd(k) = real(sum(cwd_com * grd, mask= .not. isnan(cwd_com)),r_4)
         litter_fr(k) = real(sum(litter_fr_com * grd,&
               & mask=.not.isnan(litter_fr_com)),r_4)

         cleaf_cwm(k) = real(sum(cleafcom * grd,&
               &mask=.not. isnan(cleafcom)),r_4)
         cawood_cwm(k) = real(sum(cawoodcom * grd,&
               &mask=.not. isnan(cawoodcom)),r_4)
         cfroot_cwm(k) = real(sum(cfrootcom * grd,&
               &mask=.not. isnan(cfrootcom)),r_4)

         do index = 1,6
             lnr(index,k) =  real(sum(lnr_com(index,:) * grd, &
                  & mask=.not. isnan(lnr_com(index,:))),r_4) !sdsdf
         enddo

         do index = 1,3
             storage_pool(index,k)  = real(sum(storage_pool_com(index,:) * grd,&
                  & mask=.not. isnan(storage_pool_com(index,:))),r_4) ! sdfsdf
         enddo


         ! SOIL CNP OUTPUTS
         hr_cwm(k) = real(sum(hr_com * grd, mask= .not. isnan(hr_com)),r_4)

         soil_carbon(1,k) = real(sum(soilcarbon_com(1,:) *&
                              & grd, mask= .not. isnan(soilcarbon_com(1,:))), r_4)
         soil_carbon(2,k) = real(sum(soilcarbon_com(2,:) *&
                              & grd, mask= .not. isnan(soilcarbon_com(2,:))), r_4)
         soil_carbon(3,k) = real(sum(soilcarbon_com(3,:) *&
                              & grd, mask= .not. isnan(soilcarbon_com(3,:))), r_4)
         soil_carbon(4,k) = real(sum(soilcarbon_com(4,:) *&
                              & grd, mask= .not. isnan(soilcarbon_com(4,:))), r_4)

         snr(1,k) = real(sum((snr_com(1,:) * grd), mask= .not. isnan(snr_com(1,:))), r_4)         != cwm_soil(snr(1,:), grd)
         snr(2,k) = real(sum((snr_com(2,:) * grd), mask= .not. isnan(snr_com(2,:))), r_4)         != cwm_soil(snr(2,:), grd)
         snr(3,k) = real(sum((snr_com(3,:) * grd), mask= .not. isnan(snr_com(3,:))), r_4)         != cwm_soil(snr(3,:), grd)
         snr(4,k) = real(sum((snr_com(4,:) * grd), mask= .not. isnan(snr_com(4,:))), r_4)         != cwm_soil(snr(4,:), grd)
         snr(5,k) = real(sum((snr_com(5,:) * grd), mask= .not. isnan(snr_com(5,:))), r_4)         != cwm_soil(snr(5,:), grd)
         snr(6,k) = real(sum((snr_com(6,:) * grd), mask= .not. isnan(snr_com(6,:))), r_4)         != cwm_soil(snr(6,:), grd)
         snr(7,k) = real(sum((snr_com(7,:) * grd), mask= .not. isnan(snr_com(7,:))), r_4)         != cwm_soil(snr(7,:), grd)
         snr(8,k) = real(sum((snr_com(8,:) * grd), mask= .not. isnan(snr_com(8,:))), r_4)         != cwm_soil(snr(8,:), grd)

         inorg_p(k) = real(sum((in_p_com * grd), mask=.not. isnan(in_p_com)), r_4)     != cwm_soil(in_p_com, grd)
         inorg_n(k) = real(sum((in_n_com * grd), mask=.not. isnan(in_n_com)), r_4)     != cwm_soil(in_n_com, grd)
         avail_p(k) = real(sum((av_p_com * grd), mask=.not. isnan(av_p_com)), r_4)     != cwm_soil(av_p_com, grd)
         sorbed_p(k) = real(sum((so_p_com * grd), mask=.not. isnan(so_p_com)), r_4)    != cwm_soil(so_p_com, grd)


         ! UPDATE COMMUNITY WEIGHTED VARIABLES for next day
         ! WATER
         t1ww = real(wsoil_cwm(k),r_4)
         t2ww = real(gsoil(k),r_4)
         t3ww = real(ssoil(k),r_4)

         do p = 1,npls
            wini(p) = t1ww
            gini(p) = t2ww
            sini(p) = t3ww
         enddo

         !CNP_SOIL
         in_p_in  = inorg_p(k)
         so_p_in  = sorbed_p(k)
         csoil_in = soil_carbon(:,k)
         snr_in   = snr(:,k)

         if(nupt(k) .lt. inorg_n(k)) then
             in_n_in  = inorg_n(k) - nupt(k)
             negative_n = 0.0
         else
            in_n_in = 0.0
            negative_n = inorg_n(k) - nupt(k)
         endif

         if(pupt(k) .lt. avail_p(k)) then
            av_p_in  = avail_p(k) - pupt(k)
            negative_p = 0.0
        else
           av_p_in = 0.0
           negative_p = avail_p(k) - nupt(k)
        endif

         if((k + run) .lt. 3000) then
            spnp=1
            if (in_n_in .gt. sp_available_n) in_n_in = sp_available_n
            if (av_p_in .gt. sp_available_p) then
               av_p_in = sp_available_p
               so_p_in = 0.0
               in_p_in = 0.0
            endif

            if (avail_p(k) .gt. sp_available_p)then
               avail_p(k) = sp_available_p
               inorg_p(k) = 0.0
               sorbed_p(k) = 0.0
            endif

            if (inorg_n(k) .gt. sp_available_n) inorg_n(k) = sp_available_n
         endif
            ! UPDATE DELTA CVEG POOLS FOR NEXT ROUND AND/OR LOOP
         ! UPDATE INOUTS
         dl_final = dl
         dr_final = dr
         dw_final = dw

         sp_in_p         = in_p_in
         sp_available_n  = in_n_in
         sp_available_p  = av_p_in
         sp_so_p         = so_p_in
         sp_csoil        = csoil_in
         sp_snr          = snr_in

         ! UPDATE CVEG POOLS FOR NEXT
         ! CLEAN NANs to prevent failure between cont_runs (pass only numbers to cff, clf, caf)
         !cleaf1_pft  = cleafcom
         !cawood1_pft = cawoodcom
         !cfroot1_pft = cfrootcom


         if(text_ts) then
            do npq = 1,npls
               if(gridocpcom(npq) .gt. 0.0) then
                  gridocpcom_int(npq) = 1
               else
                  gridocpcom_int(npq) = 0
               endif
            enddo

            ls = sum(gridocpcom_int)

            write(12345,1972)  k + run &
                 & ,cleaf_cwm(k)&
                 & ,cawood_cwm(k)&
                 & ,cfroot_cwm(k)&
                 & ,csoil_in(1)&
                 & ,csoil_in(2)&
                 & ,csoil_in(3)&
                 & ,csoil_in(4)&
                 & ,wsoil_cwm(k)&
                 & ,photo_cwm(k)&
                 & ,aresp_cwm(k)&
                 & ,npp_cwm(k)&
                 & ,hr_cwm(k)&
                 & ,rm_cwm(k)&
                 & ,rg_cwm(k)&
                 & ,wue(k)&
                 & ,cue(k)&
                 & ,rcm_cwm(k)&
                 & ,ls
         endif

1972  format (i12, 17(f15.6),i12)

         if (debug) then
            if(k .eq. 200) then
               close(1234)
               return
            endif
         endif


      enddo ! day_loop(nt1)

!!!!!!######/END MODEL INTEGRATION
      clf = cleaf1_pft
      caf = cawood1_pft
      cff = cfroot1_pft
      grid_area = gridocpcom!real(grd, kind=r_4)

      ! Final soil  POOLS


      if(text_ts) close(12345)

      ! IDENTIFYING PROCESS
      if (text_ts) then
         print *, ''
         print *, '---------',x,y
         print *, 'process number: ', getpid()
         print *, 'process number_from soil_dec: ',gpid()
         print *, 'soil_dec variables for this process: '
         print *, csoil_in(1:2), ' => litter_carbon'
         print *, csoil_in(3:4), ' => soil_carbon'
         print *, av_p_in, ' => labile_p_glob'
         print *, in_n_in, ' => mineral_n_glob'
      endif

   end subroutine caete_dyn

end module caete
