! Copyright 2017- LabTerra

!     This program is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
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

   subroutine caete_dyn(x,y,run,dt,w0,g0,s0,dcl,dca,dcf,prec,temp,p0,par,rhs&
        &,cleaf_ini,cawood_ini,cfroot_ini,emaxm,tsoil,photo_comm,aresp_comm&
        &,npp_comm,lai_comm,c_litter,c_soil,het_resp,rcm_comm,f51,runom_comm&
        &,evapm_comm,wsoil_comm,rm_comm,rg_comm,cleaf_comm,cawood_comm,cfroot_comm&
        &,grid_area,wue,cue,cdef,wfim,gfim,sfim,dl_final,dw_final,dr_final&
        &,clf,caf,cff,nitro_min,phop_lab,vcmax,specific_la,nupt,pupt&
        &,litter_l,cwd,litter_fr,lnr,storage_pool)

      use types
      use utils, only: gpid   => process_id
      use global_par
      use soil_dec, only: litc   => litter_carbon  ,&
                        & soic   => soil_carbon    ,&
                        & p_glob => available_p    ,&
                        & n_glob => available_n    ,&
                        & carb3  => carbon3

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

      ! podem deixar de ser inputs (chamada de fun ou como parametro global)
      real(r_4),dimension(npls),intent(in) :: cleaf_ini    ! Initial carbon content in leaves (kg m-2)
      real(r_4),dimension(npls),intent(in) :: cawood_ini   ! Initial carbon content in aboveground wood (kg m-2)
      real(r_4),dimension(npls),intent(in) :: cfroot_ini   ! Initial carbon content in fineroots (kg m-2)

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

      ! real(r_4),intent(in) :: n_mineral    ! Initial Mineral Content disponible (kg m-2)
      ! real(r_4),intent(in) :: p_labile     ! Initial Mineral Content disponible

      !     -----------------------------E N D-------------------------------

      !     -------------------------O U T P U T S---------------------------
      ! INTGRATED FROM COMMUNITY
      real(r_4),dimension(nt1),  intent(out) :: tsoil        ! soil temperature °C
      real(r_4),dimension(nt1),  intent(out) :: emaxm        ! Max.evapotranspiration (kg m-2 day-1)
      real(r_8),dimension(nt1),  intent(out) :: photo_comm    ! daily photosynthesis   (kgC m-2 year-1)
      real(r_8),dimension(nt1),  intent(out) :: aresp_comm    ! daily autotrophic res  (kgC m-2 year-1)
      real(r_8),dimension(nt1),  intent(out) :: npp_comm      ! daily net primary produ (kgC m-2 year-1)
      real(r_8),dimension(nt1),  intent(out) :: lai_comm      ! daily leaf area index m2 m-2
      real(r_4),dimension(2,nt1),intent(out) :: c_litter   ! daily litter carbon KgC m-2
      real(r_4),dimension(2,nt1),intent(out) :: c_soil     ! daily soil carbon KgC m-2
      real(r_4),dimension(nt1),  intent(out) :: het_resp     ! daily het resp  (kgC/m2)
      real(r_8),dimension(nt1),  intent(out) :: rcm_comm      ! Canopy resistence s m-1
      real(r_8),dimension(nt1),  intent(out) :: f51          ! Water stress modifier (dimensionless)
      real(r_8),dimension(nt1),  intent(out) :: runom_comm    ! Runoff kg m-2 day-1
      real(r_8),dimension(nt1),  intent(out) :: evapm_comm    ! Actual evapotranspiration kg m-2 day-1
      real(r_8),dimension(nt1),  intent(out) :: wsoil_comm    ! Soil moisture (kg m-2)
      real(r_8),dimension(nt1),  intent(out) :: wue,cue,cdef ! 1 - molCO2 m-2 s-1 (molH2O m-2 s-1)-1; 2 - npp/gpp; 3 - kgC m-2
      real(r_8),dimension(nt1),  intent(out) :: rm_comm       ! Plant (autotrophic) Maintenance respiration
      real(r_8),dimension(nt1),  intent(out) :: rg_comm       ! Plant (autotrophic) Growth respiration
      real(r_8),dimension(nt1),  intent(out) :: cleaf_comm    ! leaf biomass (KgC/m2)
      real(r_8),dimension(nt1),  intent(out) :: cawood_comm   ! aboveground wood biomass (KgC/m2)
      real(r_8),dimension(nt1),  intent(out) :: cfroot_comm   ! fine root biomass (KgC/m2)

      !  Under construction
      real(r_4),dimension(nt1),  intent(out) :: nitro_min   ! Nitrogen inorganic pool
      real(r_4),dimension(nt1),  intent(out) :: phop_lab    ! Phosphorus labile pool
      real(r_4),dimension(nt1),  intent(out) :: vcmax       ! mol m⁻² s⁻¹
      real(r_4),dimension(nt1),  intent(out) :: specific_la ! m²g⁻¹
      real(r_4),dimension(nt1),  intent(out) :: nupt        ! n plant uptake - CWM
      real(r_4),dimension(nt1),  intent(out) :: pupt        ! p plant uptake - CWM
      real(r_4),dimension(nt1),  intent(out) :: litter_l    ! CWM fluxes from vegetation to soil g m⁻² day⁻¹
      real(r_4),dimension(nt1),  intent(out) :: cwd         ! CWM
      real(r_4),dimension(nt1),  intent(out) :: litter_fr   ! CWM
      real(r_4),dimension(6,nt1),intent(out) :: lnr         ! Litter/Soil nutrient ratios
      real(r_4),dimension(3,nt1),intent(out) :: storage_pool! Rapid Auxiliary daily storage pool for carbon and nutrients

      ! COMMUNITY WIDE
      real(r_4),dimension(npls), intent(out) :: sfim           ! Kg m-2  final day water pools (snow; water; ice)
      real(r_4),dimension(npls), intent(out) :: wfim           ! Kg m-2
      real(r_4),dimension(npls), intent(out) :: gfim           ! Kg m-2

      real(r_8),dimension(npls),intent(out) :: grid_area      ! gridcell area fraction of pfts! ratio (0-1)

      real(r_8),dimension(npls),intent(out) :: clf            ! final carbon pools (cveg) for each pft (not scaled to cell area)
      real(r_8),dimension(npls),intent(out) :: caf            ! KgC m-2
      real(r_8),dimension(npls),intent(out) :: cff

      real(r_8),dimension(npls),intent(out) :: dl_final       ! delta(now, last day pool size)
      real(r_8),dimension(npls),intent(out) :: dr_final       ! KgC m-2
      real(r_8),dimension(npls),intent(out) :: dw_final

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
      ! outputs for budget (these variables usualy store monthly values from budget, when the model was stationary)
      ! Now they catch daily values (community weighted means) coming from budget
      real(r_4) :: epmes ! equal for all pfts - potential evapotranspiration(PET)(mm day-1)
      real(r_8),dimension(npls) :: rmes,phmes,smes,rcmes,f5mes
      real(r_8),dimension(npls) :: nppmes,laimes, armes!,_mies,csmes
      real(r_8),dimension(npls) :: emes,rmmes,rgmes,cuemes,c_defmes
      real(r_8),dimension(npls) :: cleafmes,cawoodmes,cfrootmes
      real(r_8),dimension(npls) :: gridocpmes
      real(r_8),dimension(npls) :: wuemes
      real(r_8),dimension(npls) :: grd
      real(r_4) :: pr,spre,ta,td,ipar,ru

      integer(i_4) :: p
      real(r_4) :: t1ww
      real(r_4) :: t2ww
      real(r_4) :: t3ww

      real(r_8) :: aux_var0_0x29a = nodata
      real(r_8) :: aux_var0_0x29b = nodata

      ! Next are auxiliary to tests
      integer(i_4),dimension(npls) :: gridocpmes_int
      logical(l_1),dimension(npls) :: gridocpmes_log
      integer(i_4) :: mypid_int, ls, npq
      character(len=29) :: filename
      character(len=7) :: mypid
      character(len=3) :: str_x, str_y
      character(len=9), parameter :: fname = 'pools_din'
      character(len=4), parameter :: fext = '.txt'
      real(r_4),dimension(8) :: soil_nr_out

      !  END OF DECLARATION


      ! ### MODEL INITIALIZATION ---------------------
      ! CALL here to initialize nutrient pools in soil

      ! store initial soil nutrients contents

      aux_var0_0x29a = n_glob
      aux_var0_0x29b = p_glob

      print*, n_glob, p_glob


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

      ! Debug mode
      if(debug) then
         open (unit=1234,file='exec_log.txt',status='unknown',&
              &form='formatted',position='append', action='write')
         write(1234,*) '--------------------------------'
         write(1234,*) '--part-1-Model initialization'
         write(1234,*) '--------------------------------'
      endif

      ! This group of variables are State Variables
      ! Inicial mass in pools ()

      ! Carbon pools initial values
      cleaf1_pft =  cleaf_ini    ! Carbon pools
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


      if(debug) then
         write(1234,*) 'cleaf1_pft - ', cleaf1_pft
         write(1234,*) 'cawood1_pft - ', cawood1_pft
         write(1234,*) 'cfroot1_pft - ', cfroot1_pft
         write(1234,*) 'dl - ', dl
         write(1234,*) 'dr - ', dr
         write(1234,*) 'dw - ', dw
         write(1234,*) '--------------------------------'
         write(1234,*) '--part-1-END'
         write(1234,*) '--------------------------------'
      endif

      !82 columns----------------------------------------------------------------------

      ! Initialize variables - Internal Variables

      cleaf_comm  = 0.0    ! leaf biomass (KgC/m2)
      cawood_comm = 0.0    ! aboveground wood biomass (KgC/m2)
      cfroot_comm = 0.0    ! fine root biomass (KgC/m2)
      grid_area  = 0.0     ! gridcell area fraction of pfts(ratio)

      emaxm     =  0.0     ! Maximum evapotranspiration
      photo_comm = 0.0     ! daily photosynthesis (kgC/m2/y)
      aresp_comm = 0.0     ! daily autotrophic respiration (kgC/m2/y)
      npp_comm   = 0.0     ! daily net primary productivity (sum of PLSs) (kgC/m2/y)
      lai_comm   = 0.0     ! daily leaf area index m2 m-2n
      ! c_litter  = 0.0    ! daily litter carbon g m-2
      ! c_soil    = 0.0    ! daily soil carbon g m-2
      het_resp  = 0.0      ! daily heterotrophic respiration (kgC/m2/y)
      rcm_comm   = 0.0     ! molH2O m-2 s-1
      f51       = 0.0      ! dimensionless
      gsoil     = 0.0      ! Soil ice mm
      ssoil     = 0.0      ! Soil snow mm
      runom_comm = 0.0     ! Runoff mm
      evapm_comm = 0.0     ! Actual evapotranspiration mm/day
      wsoil_comm = 0.0     ! Soil moisture (mm)
      rm_comm    = 0.0     ! Maintenance respiration
      rg_comm    = 0.0     ! Growth respiration
      wue       = 0.0      ! Water use efficiency (Medlyn et al. 2011 Reconciling...)
      cue       = 0.0      ! Carbon Use efficiency (NPP/GPP)

      !  ### End model initialization -------------------------------------------------

      if(debug) then
         write(1234,*) '--------------------------------'
         write(1234,*) '--part-2-daily loop'
         write(1234,*) '--------------------------------'
      endif

      !     ======================
      !     START TIME INTEGRATION
      !     ======================
      do k = 1,nt1 ! SUBSTITUTE nt1 for ndays
         if(debug) write(1234,*) '-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-LOOP: ', k

         if (k .gt. 1) then
            tsoil(k) = soil_temp(tsoil(k-1), temp(k) - 273.15)
            if(debug) write(1234,*) 'tsoil',tsoil(k), 'loop:', k
         else
            call soil_temp_sub(temp(1:1095) - 273.15,tsoil(k))
            !tsoil(k) = t0
            if(debug) write(1234,*) 'tsoil',tsoil(k), 'loop:', k
         endif

         ! if(k .eq. 1) then
         !    c_litter(:,k)  =litc
         !    c_soil(:,k)    = soic
         ! ! else
         ! !    c_litter(:,k) = c_litter(:,k-1)
         ! !    c_soil(:,k) = c_soil(:,k-1)
         ! endif
         ! ! if(k .eq. 1) then
         ! !    ! if(k .lt. 1e4) then
         ! !        n_mineral = 1
         ! !        p_labile = 0.5
         ! !    ! endif
         ! ! endif

         ! Converting units
         td = tsoil(k)
         spre = p0(k) * 0.01 ! transforamando de Pascal pra mbar (hPa)
         ta = temp(k) - 273.15 ! K to °C
         pr = prec(k) * 86400!2.62974e+06 ! kg m-2 s-1 to mm month-1
         ipar = (0.5 * par(k)) / 2.18e5 ! W m-2 to mol m-2 s-1 ! 0.5 converts RSDS to PAR
         ru = rhs(k) / 100.0
         if(debug) write(1234,*) 'ps              ta                  pr                 ipar                ru      '
         if(debug) write(1234,*) spre, ta, pr, ipar, ru

         ! Daily budget
         ! ====================
         ! Variables that I will fill with data in call_daily_budget
         epmes = 0.0
         wfim = 0.0
         gfim = 0.0
         sfim = 0.0
         smes = 0.0
         rmes = 0.0
         emes = 0.0
         phmes = 0.0
         armes = 0.0
         nppmes = 0.0
         laimes = 0.0
         f5mes = 0.0
         rcmes = 0.0
         rmmes = 0.0
         rgmes = 0.0
         cleafmes = 0.0
         cawoodmes = 0.0
         cfrootmes = 0.0
         gridocpmes = 0.0
         wuemes = 0.0
         cuemes = 0.0
         c_defmes = 0.0


         call daily_budget(dt, wini, gini,sini,td,ta,pr,spre,ipar,ru,n_glob,p_glob&
              &,storage_pool_com,cleaf1_pft,cawood1_pft,cfroot1_pft&
              &,dl,dw,dr,wfim,gfim,sfim,smes,rmes,emes,epmes,phmes,armes,nppmes&
              &,laimes,rcmes,f5mes,rmmes,rgmes,cleafmes,cawoodmes&
              &,cfrootmes,gridocpmes,wuemes,cuemes,c_defmes,vcmax_com&
              &,specific_la_com,nupt_com,pupt_com,litter_l_com,cwd_com&
              &,litter_fr_com,lnr_com)

         ! SAVE DAILY VALUES
         !82 columns---------------------------------------------------------------
         grd = gridocpmes
         emaxm(k) = epmes
         ! SOIL WATER
         gsoil(k)     = real(sum(gfim * grd, mask= .not. isnan(gfim)),r_4)
         ssoil(k)     = real(sum(sfim * grd, mask= .not. isnan(sfim)),r_4)
         wsoil_comm(k) = real(sum(wfim * grd, mask= .not. isnan(wfim)),r_4)
         snowm(k)     = real(sum(smes * grd, mask= .not. isnan(smes)),r_4)
         runom_comm(k) = real(sum(rmes * grd, mask= .not. isnan(rmes)),r_4)
         evapm_comm(k) = real(sum(emes * grd, mask= .not. isnan(emes)),r_4)
         f51(k)       = real(sum(f5mes* grd, mask= .not. isnan(f5mes)),r_4)
         ! PRODUCTIVITY
         rcm_comm(k)   = real(sum(rcmes * grd, mask= .not. isnan(rcmes)),r_4)
         lai_comm(k)   = real(sum(laimes* grd, mask= .not. isnan(laimes)),r_4)
         photo_comm(k) = real(sum(phmes * grd, mask= .not. isnan(phmes)),r_4)
         aresp_comm(k) = real(sum(armes * grd, mask= .not. isnan(armes)),r_4)
         npp_comm(k)   = real(sum(nppmes* grd, mask= .not. isnan(nppmes)),r_4)
         rm_comm(k)    = real(sum(rmmes * grd, mask= .not. isnan(rmmes)),r_4)
         rg_comm(k)    = real(sum(rgmes * grd, mask= .not. isnan(rgmes)),r_4)
         wue(k)       = real(sum(wuemes* grd, mask= .not. isnan(wuemes)),r_4)
         cue(k)       = real(sum(cuemes* grd, mask= .not. isnan(cuemes)),r_4)
         cdef(k)      = real(sum(c_defmes* grd, mask= .not. isnan(c_defmes)),r_4)

         ! Physiological traits
         vcmax(k)  = real(sum(vcmax_com* grd, mask= .not. isnan(vcmax_com)),r_4)
         specific_la(k) = real(sum(specific_la_com  * grd,&
              & mask=.not. isnan(specific_la_com)),r_4)

         ! nutrient uptake
         nupt(k) = real(sum(nupt_com * grd, mask= .not. isnan(nupt_com)),r_4)
         pupt(k) = real(sum(pupt_com * grd, mask= .not. isnan(pupt_com)),r_4)

         litter_l(k) = real(sum(litter_l_com * grd,&
              & mask=.not. isnan(litter_l_com)),r_4)
         cwd(k) = real(sum(cwd_com * grd, mask= .not. isnan(cwd_com)),r_4)
         litter_fr(k) = real(sum(litter_fr_com * grd,&
              & mask=.not.isnan(litter_fr_com)),r_4)

         cleaf_comm(k) = real(sum(cleafmes * grd,&
              &mask=.not. isnan(cleafmes)),r_4)
         cawood_comm(k) = real(sum(cawoodmes * grd,&
              &mask=.not. isnan(cawoodmes)),r_4)
         cfroot_comm(k) = real(sum(cfrootmes * grd,&
              &mask=.not. isnan(cfrootmes)),r_4)

         do index = 1,6
            lnr(index,k) =  real(sum(lnr_com(index,:) * grd, &
                 & mask=.not. isnan(lnr_com(index,:))),r_4) !sdsdf
         enddo

         do index = 1,3
            storage_pool(index,k)  = real(sum(storage_pool_com(index,:) * grd,&
                 & mask=.not. isnan(storage_pool_com(index,:))),r_4) ! sdfsdf
         enddo


         ! UPDATE WATER POOLS
         ! All PLSs receives the same value, the community weighted mean
         ! This happens to the water and mineral pools, inasmuch as these
         ! pools are shared by the community.

         t1ww = real(wsoil_comm(k),r_4)
         t2ww = real(gsoil(k),r_4)
         t3ww = real(ssoil(k),r_4)

         do p = 1,npls
            wini(p) = t1ww
            gini(p) = t2ww
            sini(p) = t3ww
         enddo

         ! UPDATE Soil Pools

         ! if(run == 0) then
         !    call carbon2(td,f51(k),evapm_comm(k),lai_comm(k),&
         !                & litc, soic)
         ! endif



         call carb3(td, (t1ww / wmax),litter_l(k),cwd(k),litter_fr(k),lnr(:,k)&
              &, litc,soic,c_litter(:,k),c_soil(:,k)&
              &,soil_nr_out, het_resp(k))

         litc = c_litter(:,k)
         soic = c_soil(:,k)

         nitro_min(k) = real(n_glob,r_4) ! - (nupt(k) * 1e-3)
         phop_lab(k) = real(p_glob,r_4) ! - (pupt(k) * 1e-3)


        ! UPDATE MINERAL POOLS
        !  n_glob = real(nitro_min(k),r_8)
        !  p_glob = real(phop_lab(k), r_8)

         ! UPDATE DELTA CVEG POOLS FOR NEXT ROUND AND/OR LOOP
         ! UPDATE INOUTS
         dl_final = dl
         dr_final = dr
         dw_final = dw

         ! UPDATE CVEG POOLS FOR NEXT LOOP
         ! CLEAN NANs to prevent failure between cont_runs (pass only numbers to cff, clf, caf)
         cleaf1_pft  = cleafmes
         cawood1_pft = cawoodmes
         cfroot1_pft = cfrootmes


        ! gridocpmes_log = gridocpmes .gt. 0.0
        ! gridocpmes_int = gridocpmes_log

         if(text_ts) then
            do npq = 1,npls
               gridocpmes_log(npq) = (gridocpmes(npq) .gt. 0.0)
            enddo

            do npq = 1, npls
               gridocpmes_int(npq) = gridocpmes_log(npq)
            enddo

            ls = sum(gridocpmes_int)

            write(12345,1972)  k + run &
                 & ,cleaf_comm(k)&
                 & ,cawood_comm(k)&
                 & ,cfroot_comm(k)&
                 & ,c_litter(1,k)&
                 & ,c_litter(2,k)&
                 & ,c_soil(1,k)&
                 & ,c_soil(2,k)&
                 & ,c_soil(3,k)&
                 & ,wsoil_comm(k)&
                 & ,photo_comm(k)&
                 & ,aresp_comm(k)&
                 & ,npp_comm(k)&
                 & ,het_resp(k)&
                 & ,rm_comm(k)&
                 & ,rg_comm(k)&
                 & ,wue(k)&
                 & ,cue(k)&
                 & ,rcm_comm(k)&
                 & ,ls
         endif

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
      grid_area = gridocpmes

      if(text_ts) close(12345)

      !  ! FINAL details
      !  ! Fill no_data in the tail of thr array (nday+1 : nt1)
      ! if(ndays < nt1) then
      !    call fill_no_data_4b(tsoil, ndays)
      !    call fill_no_data_4b(emaxm, ndays)
      !    call fill_no_data_8b(photo_comm, ndays)
      !    call fill_no_data_8b(aresp_comm, ndays)
      !    call fill_no_data_8b(npp_comm, ndays)
      !    call fill_no_data_8b(lai_comm, ndays)

      !    !# loop
      !    ! c_litter has 2 pools
      !    ! c_soil has 3 pools but the 3rd one is empty (full of zeroes)

      !    do nindex = 1, 3
      !       if(nindex < 3) call fill_no_data_4b(c_litter(nindex, :), ndays)
      !       call fill_no_data_4b(c_soil(nindex, :), ndays)
      !    enddo

      !    call fill_no_data_4b(het_resp, ndays)
      !    call fill_no_data_8b(rcm_comm, ndays)
      !    call fill_no_data_8b(f51, ndays)
      !    call fill_no_data_8b(runom_comm, ndays)
      !    call fill_no_data_8b(evapm_comm, ndays)
      !    call fill_no_data_8b(wsoil_comm, ndays)
      !    call fill_no_data_8b(wue, ndays)
      !    call fill_no_data_8b(rm_comm, ndays)
      !    call fill_no_data_8b(rg_comm, ndays)
      !    call fill_no_data_8b(cleaf_comm, ndays)
      !    call fill_no_data_8b(cawood_comm, ndays)
      !    call fill_no_data_8b(cfroot_comm, ndays)
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
         print *, litc, ' => litter_carbon'
         print *, soic, ' => soil_carbon'
         print *, p_glob, ' => labile_p_glob'
         print *, n_glob, ' => mineral_n_glob'
      endif

1972  format (i12, 18(f15.6),i12)

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


   end subroutine caete_dyn

end module caete
