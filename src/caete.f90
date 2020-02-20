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

   subroutine caete_dyn(x, y, run, dt, w0, g0, s0, csoil_init, snr_init, in_p_init, in_n_init,&
        & av_p_init, so_p_init, dcl, dca, dcf, prec, temp, p0, par, rhs,&
        & cleaf_ini, cawood_ini, cfroot_ini, emaxm, tsoil, photo_cwm, aresp_cwm, npp_cwm, lai_cwm,&
        & hr_cwm, rcm_cwm, f51, runom_cwm, evapm_cwm, wsoil_cwm, rm_cwm, rg_cwm, cleaf_cwm,&
        & cawood_cwm, cfroot_cwm, grid_area, wue, cue, cdef, wfim, gfim, sfim, dl_final, dw_final,&
        & dr_final, clf, caf, cff, vcmax, specific_la, nupt, pupt,litter_l, cwd, litter_fr,&
        & lnr, storage_pool, soil_carbon, snr, avail_p, inorg_n, inorg_p, sorbed_p)

      use types
      use utils, only: gpid   => process_id
      use global_par
      use water, only : soil_temp, soil_temp_sub
      use budget, only : daily_budget

      !     --------------------------I N P U T S----------------------------

      !boundary conditions - initial state
      integer(i_4), intent(in) :: x, y, run  ! gridcell ID
      ! integer(i_4), intent(in) :: ndays      ! Number of simulated days
      real(r_4),dimension(ntraits, npls),intent(in) :: dt ! trait table
      real(r_4),dimension(npls),intent(in) :: w0 ! Initial soil water kg m-2
      real(r_4),dimension(npls),intent(in) :: g0 ! Initial soil ice
      real(r_4),dimension(npls),intent(in) :: s0 ! Initial soil snow

      real(r_4),dimension(4),intent(in) :: csoil_init
      real(r_4),dimension(8),intent(in) :: snr_init
      real(r_4),intent(in) :: in_p_init
      real(r_4),intent(in) :: in_n_init
      real(r_4),intent(in) :: av_p_init
      real(r_4),intent(in) :: so_p_init

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

      real(r_4),dimension(npls),intent(in) :: cleaf_ini    ! Initial carbon content in leaves (kg m-2)
      real(r_4),dimension(npls),intent(in) :: cawood_ini   ! Initial carbon content in aboveground wood (kg m-2)
      real(r_4),dimension(npls),intent(in) :: cfroot_ini   ! Initial carbon content in fineroots (kg m-2)

      !     -----------------------------E N D-------------------------------

      !     -------------------------O U T P U T S---------------------------

      real(r_4),dimension(nt1),   intent(out) :: emaxm        ! 1   ! Max.evapotranspiration (kg m-2 day-1)
      real(r_4),dimension(nt1),   intent(out) :: tsoil        ! 2   ! soil temperature °C
      real(r_8),dimension(nt1),   intent(out) :: photo_cwm    ! 3   ! daily photosynthesis   (kgC m-2 year-1)
      real(r_8),dimension(nt1),   intent(out) :: aresp_cwm    ! 4   ! daily autotrophic res  (kgC m-2 year-1)
      real(r_8),dimension(nt1),   intent(out) :: npp_cwm      ! 5   ! daily net primary produ (kgC m-2 year-1)
      real(r_8),dimension(nt1),   intent(out) :: lai_cwm      ! 6   ! daily leaf area index m2 m-2
      real(r_4),dimension(nt1),   intent(out) :: hr_cwm       ! 7   ! daily het resp  (kgC/m2)
      real(r_8),dimension(nt1),   intent(out) :: rcm_cwm      ! 8  ! leaf resistence s m-1
      real(r_8),dimension(nt1),   intent(out) :: f51          ! 9  ! Water stress modifier (dimensionless)
      real(r_8),dimension(nt1),   intent(out) :: runom_cwm    ! 10  ! Runoff kg m-2 day-1
      real(r_8),dimension(nt1),   intent(out) :: evapm_cwm    ! 11  ! Actual evapotranspiration kg m-2 day-1
      real(r_8),dimension(nt1),   intent(out) :: wsoil_cwm    ! 12  ! Soil moisture (kg m-2)
      real(r_8),dimension(nt1),   intent(out) :: rm_cwm       ! 13  ! Plant (autotrophic) Maintenance respiration
      real(r_8),dimension(nt1),   intent(out) :: rg_cwm       ! 14  ! Plant (autotrophic) Growth respiration
      real(r_8),dimension(nt1),   intent(out) :: cleaf_cwm    ! 15  ! leaf biomass (KgC/m2)
      real(r_8),dimension(nt1),   intent(out) :: cawood_cwm   ! 16  ! aboveground wood biomass (KgC/m2)
      real(r_8),dimension(nt1),   intent(out) :: cfroot_cwm   ! 17  ! fine root biomass (KgC/m2)
      real(r_8),dimension(npls),  intent(out) :: grid_area    ! 18  ! gridcell area fraction of pfts! ratio (0-1)
      real(r_8),dimension(nt1),   intent(out) :: wue          ! 19  ! 1 - molCO2 m-2 s-1 (molH2O m-2 s-1)-1;
      real(r_8),dimension(nt1),   intent(out) :: cue          ! 20  ! 2 - npp/gpp;
      real(r_8),dimension(nt1),   intent(out) :: cdef         ! 21  ! 3 - kgC m-2
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

      real(r_4),dimension(6,nt1), intent(out) :: lnr          ! 38 ! Litter nutrient ratios
      real(r_4),dimension(3,nt1), intent(out) :: storage_pool ! 39 ! Rapid Auxiliary daily storage pool for carbon and nutrients g m⁻²

      real(r_4),dimension(4,nt1), intent(out) :: soil_carbon  ! 40 ! CSOIL_POOL carbon gC m-2
      real(r_4),dimension(8,nt1), intent(out) :: snr          ! 41 ! Soil nutrient ratios g/g
      real(r_4),dimension(nt1),   intent(out) :: avail_p      ! 42
      real(r_4),dimension(nt1),   intent(out) :: inorg_n      ! 43
      real(r_4),dimension(nt1),   intent(out) :: inorg_p      ! 44
      real(r_4),dimension(nt1),   intent(out) :: sorbed_p     ! 45

      !  c     --------------------------------E N D----------------------------

      !  c     ------------------------- internal variables---------------------
      integer(i_4) :: k, index, nindex
      real(r_8),dimension(nt1) :: gsoil    !Soil ice kg m-2
      real(r_8),dimension(nt1) :: ssoil    !Soil snow kg m-2
      real(r_8),dimension(nt1) :: snowm    !Snowmelt kg m-2

      real(r_8),dimension(3,npls) :: storage_pool_com
      real(r_8),dimension(npls) :: vcmax_com
      real(r_8),dimension(npls) :: specific_la_com
      real(r_8),dimension(npls) :: nupt_com
      real(r_8),dimension(npls) :: pupt_com
      real(r_8),dimension(npls) :: litter_l_com
      real(r_8),dimension(npls) :: cwd_com
      real(r_8),dimension(npls) :: litter_fr_com
      real(r_8),dimension(6,npls) :: lnr_com

      real(r_4),dimension(npls) :: sini  !kg m-2
      real(r_4),dimension(npls) :: wini  !kg m-2
      real(r_4),dimension(npls) :: gini  !kg m-2

      real(r_8),dimension(npls) :: dl  ! delta(now, last day pool size)
      real(r_8),dimension(npls) :: dr  !kg m-2
      real(r_8),dimension(npls) :: dw

      real(r_8),dimension(npls) :: cleaf1_pft   !kg m-2
      real(r_8),dimension(npls) :: cawood1_pft  !kg m-2
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
      real(r_8),dimension(npls) :: hr_com
      real(r_8),dimension(npls) :: wue_com

      real(r_4),dimension(4,npls) :: soilcarbon_com
      real(r_4),dimension(8,npls) :: snr_com
      real(r_4),dimension(npls)   :: in_p_com
      real(r_4),dimension(npls)   :: in_n_com
      real(r_4),dimension(npls)   :: av_p_com
      real(r_4),dimension(npls)   :: so_p_com
      real(r_4) :: csoil_in, snr_in, in_p_in, in_n_in, av_p_in, so_p_in ! Auxiliary variables
      real(r_8),dimension(npls) :: grd
      real(r_4) :: pr,spre,ta,td,ipar,ru
      integer(i_4) :: p
      real(r_4) :: t1ww
      real(r_4) :: t2ww
      real(r_4) :: t3ww

      real(r_4) :: aux1, aux2
      ! Next are auxiliary to tests
      integer(i_4),dimension(npls) :: gridocpcom_int
      logical(l_1),dimension(npls) :: gridocpcom_log
      integer(i_4) :: mypid_int, ls, npq
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
      in_p_in  = in_p_init
      in_n_in  = in_n_init
      av_p_in  = av_p_init
      so_p_in  = so_p_init
      csoil_in = csoil_init
      snr_in   = snr_init

      !82 columns-----------------------------------------------------------------

      ! Initialize variables
      cleaf_cwm          = 0.0     ! leaf biomass (KgC/m2)
      cawood_cwm         = 0.0     ! aboveground wood biomass (KgC/m2)
      cfroot_cwm         = 0.0     ! fine root biomass (KgC/m2)
      grid_area          = 0.0     ! occupation coefficient
      emaxm              = 0.0     ! Maximum evapotranspiration
      photo_cwm          = 0.0     ! daily photosynthesis (kgC/m2/y)
      aresp_cwm          = 0.0     ! daily autotrophic respiration (kgC/m2/y)
      npp_cwm            = 0.0     ! daily net primary productivity (sum of PLSs) (kgC/m2/y)
      lai_cwm            = 0.0     ! daily leaf area index m2 m-2n
      csoil              = 0.0     ! Carbon in soil g C m-2 (l1, l2, s1, s2)
      hr_cwm             = 0.0     ! daily heterotrophic respiration (kgC/m2/y)
      rcm_cwm            = 0.0     ! molH2O m-2 s-1
      f51                = 0.0     ! dimensionless
      gsoil              = 0.0     ! Soil ice mm
      ssoil              = 0.0     ! Soil snow mm
      runom_cwm          = 0.0     ! Runoff mm
      evapm_cwm          = 0.0     ! Actual evapotranspiration mm/day
      wsoil_cwm          = 0.0     ! Soil moisture (mm)
      rm_cwm             = 0.0     ! Maintenance respiration
      rg_cwm             = 0.0     ! Growth respiration
      wue                = 0.0     ! Water use efficiency (Medlyn et al. 2011 Reconciling...)
      cue                = 0.0     ! Carbon Use efficiency (NPP/GPP)
      soilcarbon_com     = 0.0
      snr_com            = 0.0
      in_p_com           = 0.0
      in_n_com           = 0.0
      av_p_com           = 0.0
      so_p_com           = 0.0

      !  ### End model initialization -------------------------------------------------


      !     ======================
      !     START TIME INTEGRATION
      !     ======================
      do k = 1,nt1
         if (k .gt. 1) then
            tsoil(k) = soil_temp(tsoil(k-1), temp(k) - 273.15)
            if(debug) write(1234,*) 'tsoil',tsoil(k), 'loop:', k
         else
            call soil_temp_sub(temp(1:1095) - 273.15,tsoil(k))
            !tsoil(k) = t0
            if(debug) write(1234,*) 'tsoil',tsoil(k), 'loop:', k
         endif

         ! Setting input variables and Converting units
         td   = tsoil(k)
         spre = p0(k) * 0.01 ! transforamando de Pascal pra mbar (hPa)
         ta   = temp(k) - 273.15 ! K to °C
         pr   = prec(k) * 86400!2.62974e+06 ! kg m-2 s-1 to mm month-1
         ipar = (0.5 * par(k)) / 2.18e5 ! W m-2 to mol m-2 s-1 ! 0.5 converts RSDS to PAR
         ru   = rhs(k) / 100.0

         ! Daily budget
         ! ====================
         ! Variables that I will fill with data in call daily_budget
         ep_daily                = 0.0
         wfim                    = 0.0
         gfim                    = 0.0
         sfim                    = 0.0
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
         cleafcom                = 0.0
         cawoodcom               = 0.0
         cfrootcom               = 0.0
         gridocpcom              = 0.0
         wue_com                 = 0.0
         cue_com                 = 0.0
         c_def_com               = 0.0
         hr_com                  = 0.0


         ! ! AVAILABLE NUTRIENTS:
         ! if (run .lt. 5000) then
         !    av_n = available_n_init
         !    av_p = available_p_init
         ! else
         !    av_n = available_n
         !    av_p = available_p
         ! endif

         ! subroutine daily_budget(dt, w1, g1, s1, ts, temp, prec, p0, ipar, rh,&
         !    & inorg_p, inorg_n, avail_p, sorb_p, sto_budg,&
         !    & cl1_pft, ca1_pft, cf1_pft, dleaf, dwood, droot,&
         !    & csoil, snr_com, w2, g2, s2, smavg, ruavg, evavg, epavg,&
         !    & phavg, aravg, nppavg, laiavg, rcavg, f5avg,&
         !    & rmavg, rgavg, cleafavg_pft, cawoodavg_pft, cfrootavg_pft,&
         !    & ocpavg, wueavg, cueavg, c_defavg, vcmax, specific_la, soilc,&
         !    & inorganic_p, inorganic_n, available_p, sorbed_p,&
         !    & nupt, pupt, litter_l, cwd, litter_fr, het_resp, lnr, snr)


         call daily_budget(dt, wini, gini, sini, td, ta, pr, spre, ipar, ru&
              &, in_p_in, in_n_in, av_p_in, so_p_in, storage_pool_com, cleaf1_pft&
              &, cawood1_pft, cfroot1_pft, dl, dw, dr, csoil_in, snr_in, wfim, gfim, sfim&
              &, snowmelt_com, runoff_com, e_com, ep_daily, gpp_com, ar_com, npp_com&
              &, lai_com, canopy_res_com, f5_com, rm_com, rg_com, cleafcom, cawoodcom&
              &, cfrootcom, gridocpcom, wue_com, cue_com, c_def_com, vcmax_com&
              &, specific_la_com, soilcarbon_com, in_p_com, in_n_com, av_p_com, so_p_com&
              &, nupt_com, pupt_com, litter_l_com, cwd_com&
              &, litter_fr_com, hr_com, lnr_com, snr_com)

              ! TO DO
              ! REZA A LENDA QUE TUDO QUE FOR MES VAI VIRAR DAILY
              ! CRIAR AS NOVAS VARIAVEIS NESTE ESCOPO
              ! APLICAR A CWM e limpar a porra toda

         !82 columns-------------------------------------------------------------
         grd            = gridocpcom
         emaxm(k)       = ep_daily
         gsoil(k)       = cwm(real(gfim, r_8), grd)
         ssoil(k)       = cwm(real(sfim, r_8), grd)
         wsoil_cwm(k)   = cwm(real(wfim, r_8), grd)
         snowm(k)       = cwm(snowmelt_com, grd)
         runom_cwm(k)   = cwm(runoff_com, grd)
         evapm_cwm(k)   = cwm(e_com, grd)
         f51(k)         = cwm(f5_com, grd)
         rcm_cwm(k)     = cwm(canopy_res_com, grd)
         lai_cwm(k)     = cwm(lai_com, grd)
         photo_cwm(k)   = cwm(gpp_com, grd)
         aresp_cwm(k)   = cwm(ar_com, grd)
         npp_cwm(k)     = cwm(npp_com, grd)
         rm_cwm(k)      = cwm(rm_com, grd)
         rg_cwm(k)      = cwm(rg_com, grd)
         wue(k)         = cwm(wue_com, grd)
         cue(k)         = cwm(cue_com, grd)
         cdef(k)        = cwm(c_def_com, grd)
         vcmax(k)       = cwm(vcmax_com, grd)
         specific_la(k) = cwm(specific_la_com, grd)
         nupt(k)        = cwm(nupt_com, grd)
         pupt(k)        = cwm(pupt_com, grd)
         litter_l(k)    = cwm(litter_l_com, grd)
         cwd(k)         = cwm(cwd_com, grd)
         litter_fr(k)   = cwm(litter_fr_com, grd)
         cleaf_cwm(k)   = cwm(cleafcom, grd)
         cawood_cwm(k)  = cwm(cawoodcom, grd)
         cfroot_cwm(k)  = cwm(cfrootcom, grd)
         hr_cwm(k)      = cwm(hr_com, grd)

         ! SOIL CNP OUTPUTS
         soil_carbon(1,k) = cwm_soil(soilcarbon_com(1,:), grd)
         soil_carbon(2,k) = cwm_soil(soilcarbon_com(2,:), grd)
         soil_carbon(3,k) = cwm_soil(soilcarbon_com(3,:), grd)
         soil_carbon(4,k) = cwm_soil(soilcarbon_com(4,:), grd)

         snr(1,k)         = cwm_soil(snr(1,:), grd)
         snr(2,k)         = cwm_soil(snr(2,:), grd)
         snr(3,k)         = cwm_soil(snr(3,:), grd)
         snr(4,k)         = cwm_soil(snr(4,:), grd)
         snr(5,k)         = cwm_soil(snr(5,:), grd)
         snr(6,k)         = cwm_soil(snr(6,:), grd)
         snr(7,k)         = cwm_soil(snr(7,:), grd)
         snr(8,k)         = cwm_soil(snr(8,:), grd)

         inorg_p(k)     = cwm_soil(in_p_com, grd)
         inorg_n(k)     = cwm_soil(in_n_com, grd)
         avail_p(k)     = cwm_soil(av_p_com, grd)
         sorbed_p(k)    = cwm_soil(so_p_com, grd)


         do index = 1,6
            lnr(index,k) = cwm(lnr_com(index,:), grd)
         enddo

         do index = 1,3
            storage_pool(index,k)  = cwm(storage_pool_com(index,:), grd)
         enddo

         ! UPDATE COMMUNITY WEIGHTED VARIABLES
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
         in_n_in  = inorg_n(k)
         av_p_in  = avail_p(k)
         so_p_in  = sorbed_p(k)
         csoil_in = soil_carbon(:,k)
         snr_in   = snr(:,k)


         ! UPDATE DELTA CVEG POOLS FOR NEXT ROUND AND/OR LOOP
         ! UPDATE INOUTS
         dl_final = dl
         dr_final = dr
         dw_final = dw

         ! UPDATE CVEG POOLS FOR NEXT
         ! CLEAN NANs to prevent failure between cont_runs (pass only numbers to cff, clf, caf)
         cleaf1_pft  = cleafcom
         cawood1_pft = cawoodcom
         cfroot1_pft = cfrootcom

         if(text_ts) then
            do npq = 1,npls
               gridocpcom_log(npq) = (gridocpcom(npq) .gt. 0.0)
            enddo

            do npq = 1, npls
               gridocpcom_int(npq) = gridocpcom_log(npq)
            enddo

            ls = sum(gridocpcom_int)

            write(12345,1972)  k + run &
                 & ,cleaf_cwm(k)&
                 & ,cawood_cwm(k)&
                 & ,cfroot_cwm(k)&
                 & ,csoil(1,k)&
                 & ,csoil(2,k)&
                 & ,csoil(3,k)&
                 & ,csoil(4,k)&
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

      ! Change this section to prevent NANs

      clf = cleaf1_pft
      caf = cawood1_pft
      cff = cfroot1_pft
      grid_area = gridocpcom


      if(text_ts) close(12345)

      !  ! FINAL details
      !  ! Fill no_data in the tail of thr array (nday+1 : nt1)
      ! if(ndays < nt1) then
      !    call fill_no_data_4b(tsoil, ndays)
      !    call fill_no_data_4b(emaxm, ndays)
      !    call fill_no_data_8b(photo_cwm, ndays)
      !    call fill_no_data_8b(aresp_cwm, ndays)
      !    call fill_no_data_8b(npp_cwm, ndays)
      !    call fill_no_data_8b(lai_cwm, ndays)

      !    !# loop
      !    ! c_litter has 2 pools
      !    ! c_soil has 3 pools but the 3rd one is empty (full of zeroes)

      !    do nindex = 1, 3
      !       if(nindex < 3) call fill_no_data_4b(c_litter(nindex, :), ndays)
      !       call fill_no_data_4b(c_soil(nindex, :), ndays)
      !    enddo

      !    call fill_no_data_4b(hr_cwm, ndays)
      !    call fill_no_data_8b(rcm_cwm, ndays)
      !    call fill_no_data_8b(f51, ndays)
      !    call fill_no_data_8b(runom_cwm, ndays)
      !    call fill_no_data_8b(evapm_cwm, ndays)
      !    call fill_no_data_8b(wsoil_cwm, ndays)
      !    call fill_no_data_8b(wue, ndays)
      !    call fill_no_data_8b(rm_cwm, ndays)
      !    call fill_no_data_8b(rg_cwm, ndays)
      !    call fill_no_data_8b(cleaf_cwm, ndays)
      !    call fill_no_data_8b(cawood_cwm, ndays)
      !    call fill_no_data_8b(cfroot_cwm, ndays)
      !    call fill_no_data_4b(nitro_min, ndays)
      !    call fill_no_data_4b(phop_lab, ndays)
      !    call fill_no_data_4b(vcmax, ndays)
      !    call fill_no_data_4b(specific_la, ndays)
      !    call fill_no_data_4b(nupt, ndays)
      !    call fill_no_data_4b(pupt, ndays)
      !    call fill_no_data_4b(litter_l, ndays)
      !    call fill_no_data_4b(cwd, ndays)
      !    call fill_no_data_4b(litter_fr, ndays)

      !    ! Loop over sdditional dimension (lnr 6,nt1) (storage_pool 3,nt1)
      !    do nindex = 1,6
      !       if(nindex < 4) call fill_no_data_4b(storage_pool(nindex,:), ndays)
      !       call fill_no_data_4b(lnr(nindex,:), ndays)
      !    enddo
      ! endif
      ! IDENTIFYING PROCESS
      if (text_ts) then
         print *, ''
         print *, '---------',x,y
         print *, 'process number: ', getpid()
         print *, 'process number_from soil_dec: ',gpid()
         print *, 'soil_dec variables for this process: '
         print *, soc_soil_dec(1:2), ' => litter_carbon'
         print *, soc_soil_dec(3:4), ' => soil_carbon'
         print *, available_p, ' => labile_p_glob'
         print *, available_n, ' => mineral_n_glob'
      endif

   contains
         subroutine fill_no_data_4b(arr, nd)

            use types
            use global_par

            real(kind=r_4), dimension(nt1), intent(inout) :: arr
            integer(kind=i_4), intent(in) :: nd

            arr(nd + 1: nt1) = -9999.0

         end subroutine fill_no_data_4b


         subroutine fill_no_data_8b(arr, nd)

            use types
            use global_par

            real(kind=r_8), dimension(nt1), intent(inout) :: arr
            integer(kind=i_4), intent(in) :: nd

            arr(nd + 1: nt1) = -9999.0D0

         end subroutine fill_no_data_8b


         function cwm(var_arr, area_arr) result(retval)

            use types
            use global_par

            real(kind=r_8), dimension(npls), intent(in) :: var_arr, area_arr
            real(kind=r_4) :: retval
            retval = real(sum(var_arr * area_arr, mask = .not. isnan(var_arr)), r_4)

         end function cwm

         function cwm_soil(var_arr, area_arr) result(retval)

            use types
            use global_par

            real(kind=r_8), dimension(npls), intent(in) :: area_arr
            real(kind=r_4), dimension(npls), intent(in) :: var_arr
            real(kind=r_4) :: retval
            retval = sum(var_arr * real(area_arr, r_4), mask = .not. isnan(var_arr))

         end function cwm

   end subroutine caete_dyn

end module caete
