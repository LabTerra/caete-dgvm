! Copyright 2017- LabTerra

!     This program is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.

!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.

!     You should have received a copy of the GNU General Public License
!     along with this program.  If not, see <http://www.gnu.org/licenses/>.

! TODO: these changes can improve performance (memory and speed)
! CHECK the logic of THe inputs of soil water and organic matter
! INCLUDE new variables in hell
! check the logic of no_pls_run (important for performance)
! check the logic of no_pls (idem)
module budget
   implicit none
   private

   public :: daily_budget

contains

   subroutine daily_budget(dt, w1, g1, s1, ts, temp, prec, p0, ipar, rh,&
        & inorg_p, inorg_n, avail_p, sorb_p, sto_budg,&
        & cl1_pft, ca1_pft, cf1_pft, dleaf, dwood, droot,&
        & csoil_com, snr_com, w2, g2, s2, smavg, ruavg, evavg, epavg,&
        & phavg, aravg, nppavg, laiavg, rcavg, f5avg,&
        & rmavg, rgavg, cleafavg_pft, cawoodavg_pft, cfrootavg_pft,&
        & ocpavg, wueavg, cueavg, c_defavg, vcmax, specific_la, soilc,&
        & inorganic_p, inorganic_n, available_p, sorbed_p,&
        & nupt, pupt, litter_l, cwd, litter_fr, het_resp, lnr, snr)

      use types
      use global_par
      use photo, only: pft_area_frac, allocation
      use water, only: evpot2, penman, available_energy, runoff
      use productivity
      use soil_dec, only: carb3 => carbon3
      !     ----------------------------INPUTS-------------------------------
      real(r_4),dimension(ntraits,npls),intent(in) :: dt
      real(r_4),dimension(npls),intent(in) :: w1   !Initial (previous day) soil moisture storage (mm)
      real(r_4),dimension(npls),intent(in) :: g1   !Initial soil ice storage (mm)
      real(r_4),dimension(npls),intent(in) :: s1   !Initial overland snow storage (mm)
      real(r_4),intent(in) :: ts                   ! Soil temperature (oC)
      real(r_4),intent(in) :: temp                 ! Surface air temperature (oC)
      real(r_4),intent(in) :: prec                 ! Precipitation (mm/day)
      real(r_4),intent(in) :: p0                   ! Surface pressure (mb)
      real(r_4),intent(in) :: ipar                 ! Incident photosynthetic active radiation mol Photons m-2 s-1
      real(r_4),intent(in) :: rh                   ! Relative humidity
      real(r_4),intent(in) :: inorg_p
      real(r_4),intent(in) :: inorg_n
      real(r_4),intent(in) :: avail_p
      real(r_4),intent(in) :: sorb_p
      real(r_4),dimension(4),intent(in) :: csoil_com      ! Soil carbon (gC/m2)   State Variable -> The size of the carbon pools
      real(r_4),dimension(8),intent(in)  :: snr_com

      real(r_8),dimension(3,npls),intent(inout)  :: sto_budg ! Rapid Storage Pool (C,N,P)
      real(r_8),dimension(npls),intent(inout) :: cl1_pft  ! initial BIOMASS cleaf compartment
      real(r_8),dimension(npls),intent(inout) :: ca1_pft  !                 cawood
      real(r_8),dimension(npls),intent(inout) :: cf1_pft  !                 froot

      real(r_8),dimension(npls),intent(inout) :: dleaf  ! CHANGE IN cVEG (DAILY BASIS) TO GROWTH RESP
      real(r_8),dimension(npls),intent(inout) :: dwood
      real(r_8),dimension(npls),intent(inout) :: droot

      !     ----------------------------OUTPUTS------------------------------
      real(r_4),intent(out) :: epavg                          !Maximum evapotranspiration (mm/day)
      real(r_4),intent(out),dimension(npls) :: w2             !Final (last day) soil moisture storage (mm)
      real(r_4),intent(out),dimension(npls) :: g2             !Final soil ice storage (mm)
      real(r_4),intent(out),dimension(npls) :: s2             !Final overland snow storage (mm)
      real(r_8),intent(out),dimension(npls) :: smavg          !Snowmelt Daily average (mm/day)
      real(r_8),intent(out),dimension(npls) :: ruavg          !Runoff Daily average (mm/day)
      real(r_8),intent(out),dimension(npls) :: evavg          !Actual evapotranspiration Daily average (mm/day)
      real(r_8),intent(out),dimension(npls) :: phavg          !Daily photosynthesis (Kg m-2 day-1)
      real(r_8),intent(out),dimension(npls) :: aravg          !Daily autotrophic respiration (Kg m-2 day-1)
      real(r_8),intent(out),dimension(npls) :: nppavg         !Daily NPP (average between PFTs)(Kg m-2 day-1)
      real(r_8),intent(out),dimension(npls) :: laiavg         !Daily leaf area Index m2m-2
      real(r_8),intent(out),dimension(npls) :: rcavg          !Daily canopy resistence s/m
      real(r_8),intent(out),dimension(npls) :: f5avg          !Daily canopy resistence s/m
      real(r_8),intent(out),dimension(npls) :: rmavg,rgavg    !maintenance/growth respiration (Kg m-2 day-1)
      real(r_8),intent(out),dimension(npls) :: cleafavg_pft   !Carbon in plant tissues (kg m-2)
      real(r_8),intent(out),dimension(npls) :: cawoodavg_pft  !
      real(r_8),intent(out),dimension(npls) :: cfrootavg_pft  !
      real(r_8),intent(out),dimension(npls) :: ocpavg         ! [0-1]
      real(r_8),intent(out),dimension(npls) :: wueavg         !
      real(r_8),intent(out),dimension(npls) :: cueavg         ! [0-1]
      real(r_8),intent(out),dimension(npls) :: c_defavg       ! kg(C) m-2
      real(r_8),intent(out),dimension(npls) :: vcmax          ! µmol m-2 s-1
      real(r_8),intent(out),dimension(npls) :: specific_la    ! m2 g(C)-1
      real(r_4),intent(out),dimension(4,npls) :: soilc        ! Soil carbon pools (gC m-2)
      real(r_4),intent(out),dimension(npls) :: inorganic_p
      real(r_4),intent(out),dimension(npls) :: inorganic_n
      real(r_4),intent(out),dimension(npls) :: available_p
      real(r_4),intent(out),dimension(npls) :: sorbed_p
      real(r_8),intent(out),dimension(npls) :: nupt           ! gN m-2 ! Nitrogen uptake
      real(r_8),intent(out),dimension(npls) :: pupt           ! gP m-2 ! Phosphoruns uptake
      real(r_8),intent(out),dimension(npls) :: litter_l       ! gC m-2 ! Litter from leaves
      real(r_8),intent(out),dimension(npls) :: cwd            ! gC m-2 ! coarse wood debris
      real(r_8),intent(out),dimension(npls) :: litter_fr      ! gC m-2 ! litter from fine roots
      real(r_4),intent(out),dimension(npls) :: het_resp       ! gC m-2
      ! Litter Nutrient Ratisos :: variables(6)         [(lln2c),(rln2c),(cwdn2c),(llp2c),(rlp2c),(cwdp2c)]
      real(r_8),intent(out),dimension(6,npls) :: lnr         ! g(N) g(C)-1 Litter Nutrient to C ratios (Comming from cveg pools)
      ! Soil Nutrient Ratios :: variables(8)            [(l1n2c),(l2n2c),(c1dn2c),(c2n2c),(l1p2c),(l2p2c),(c1p2c),(c2p2c)]
      real(r_4),intent(out),dimension(8,npls) :: snr         ! g(N) g(C)-1 Soil Nutrient to C ratios (IN soil C pools)
      !     -----------------------Internal Variables------------------------

      integer(i_4) :: p, index
      logical :: end_pls = .false.
      logical :: no_pls_run = .false.
      real(r_4),dimension(ntraits) :: dt1 ! Store pls attributes array (1D)
      !     RELATED WITH GRIDCELL OCUPATION

      real(r_8),dimension(npls) :: ocp_mm
      real(r_8),dimension(npls) :: ocp_coeffs !,ocp_coeffs2
      logical(l_1),dimension(npls) :: light_limitation_bool !, ocp_wood2

      !     WBM COMMUNICATION (water balance)
      real(r_4) :: psnow                !Snowfall (mm/day)
      real(r_4) :: prain                !Rainfall (mm/day)
      real(r_4) :: emax

      real(r_4),dimension(npls) :: rimelt               !Runoff due to soil ice melting
      real(r_4),dimension(npls) :: smelt                !Snowmelt (mm/day)
      real(r_4),dimension(npls) :: w                    !Daily soil moisture storage (mm)
      real(r_4),dimension(npls) :: g                    !Daily soil ice storage (mm)
      real(r_4),dimension(npls) :: s                    !Daily overland snow storage (mm)
      real(r_4),dimension(npls) :: ds
      real(r_4),dimension(npls) :: dw
      real(r_4),dimension(npls) :: roff                 !Total runoff
      real(r_4),dimension(npls) :: evap                !Actual evapotranspiration (mm/day)


      !c     Carbon Cycle
      real(r_4),dimension(npls) ::  ph             !Canopy gross photosynthesis (kgC/m2/yr)
      real(r_4),dimension(npls) ::  ar             !Autotrophic respiration (kgC/m2/yr)
      real(r_4),dimension(npls) ::  nppa           !Net primary productivity / auxiliar
      real(r_4),dimension(npls) ::  laia           !Leaf area index (m2 leaf/m2 area)
      real(r_4),dimension(npls) ::  rc2            !Canopy resistence (s/m)
      real(r_4),dimension(npls) ::  f1             !
      real(r_4),dimension(npls) ::  f5             !Photosynthesis (mol/m2/s)
      real(r_4),dimension(npls) ::  vpd            !Vapor Pressure deficit
      real(r_4),dimension(npls) ::  rm             ! maintenance & growth a.resp
      real(r_4),dimension(npls) ::  rg
      real(r_4),dimension(npls) ::  wue, cue, c_def
      real(r_8),dimension(npls) ::  cl1,cf1,ca1 ! carbon pre-allocation
      real(r_8),dimension(npls) ::  cl2,cf2,ca2 ! carbon pos-allocation
      real(r_8),dimension(3,npls) :: day_storage   ! g m-2
      real(r_8),dimension(npls) :: n_uptake           ! g m-2
      real(r_8),dimension(npls) :: p_uptake           ! g m-2
      ! CNP CYCLE VARIABLES
      real(r_4),dimension(2, npls) :: litter_carbon_bdg, soil_carbon_bdg, clitter, csoil
      real(r_4),dimension(8, npls) :: snr_internal = 0.0
      real(r_4),dimension(8, npls) :: snr_aux = 0.0
      real(r_4),dimension(npls) :: in_p, in_n, av_p, so_p

      real(r_8) :: litter_fr_aux, litter_l_aux, cwd_aux

      !     Precipitation
      !     =============
      psnow = 0.0
      prain = 0.0
      if (temp.lt.tsnow) then
         psnow = prec !* 0.0328549 ! converts mm/month in mm/day
      else
         prain = prec !* 0.0328549
      endif

      !     Initialization
      !     --------------
      !  plss vectors (outputs)
      epavg                 = 0.0
      smavg                 = 0.0
      ruavg                 = 0.0
      evavg                 = 0.0
      rcavg                 = 0.0
      laiavg                = 0.0
      phavg                 = 0.0
      aravg                 = 0.0
      nppavg                = 0.0
      rmavg                 = 0.0
      rgavg                 = 0.0
      ocpavg                = 0.0
      wueavg                = 0.0
      cueavg                = 0.0
      ocp_mm                = 0.0
      emax                  = 0.0
      nupt                  = 0.0
      pupt                  = 0.0
      litter_l              = 0.0
      cwd                   = 0.0
      litter_fr             = 0.0
      het_resp              = 0.0
      nppa                  = 0.0
      ph                    = 0.0
      ar                    = 0.0
      laia                  = 0.0
      f5                    = 0.0
      f1                    = 0.0
      vpd                   = 0.0
      rc2                   = 0.0
      rm                    = 0.0
      rg                    = 0.0
      wue                   = 0.0
      cue                   = 0.0
      rc2                   = 0.0
      f1                    = 0.0
      f5                    = 0.0
      ocp_coeffs            = 0.0D0
      light_limitation_bool = .false.

      ! Persistent pools
      w                     = w1     ! hidrological pools state vars
      g                     = g1
      s                     = s1
      cl1                   = cl1_pft ! daily initial carbonVEG pools
      ca1                   = ca1_pft
      cf1                   = cf1_pft

      do p = 1, npls
         clitter(1,p)    = csoil_com(1)
         clitter(2,p)    = csoil_com(2)
         csoil(1,p)      = csoil_com(3)
         csoil(2,p)      = csoil_com(4)
         in_p(p)                   = inorg_p
         in_n(p)                   = inorg_n
         av_p(p)                   = avail_p
         so_p(p)                   = sorb_p
         do index = 1, 8
            snr_internal(index,p)  = snr_com(index)
         enddo
      enddo
      !     Grid cell area fraction (%) ocp_coeffs(pft(1), pft(2), ...,pft(p))
      !     =================================================================
      call pft_area_frac(cl1, cf1, ca1, ocp_coeffs, light_limitation_bool) ! def in funcs.f90

      !     Maximum evapotranspiration   (emax)
      !     =================================
      emax = evpot2(p0,temp,rh,available_energy(temp))

      !     Productivity & Growth (ph, ALLOCATION, aresp, vpd, rc2 & etc.) for each PLS
      !     =================================
      do p = 1,npls

         dt1 = dt(:,p) ! Pick up the pls functional attributes list of PLS npls
         litter_fr_aux = 0.0D0
         litter_l_aux = 0.0D0
         cwd_aux = 0.0D0

         end_pls = .false.

         call prod(dt1, light_limitation_bool(p), temp, p0, w(p), ipar, rh, emax, cl1(p)&
              &, ca1(p), cf1(p), dleaf(p), dwood(p), droot(p), sto_budg(:,p), ph(p), ar(p)&
              &, nppa(p),laia(p), f5(p), vpd(p), rm(p), rg(p), rc2(p), wue(p), c_def(p)&
              &, vcmax(p), specific_la(p), day_storage(:,p))

         sto_budg(:,p) = day_storage(:,p)

         !     Carbon/Nitrogen/Phosphorus allocation/deallocation
         !     =====================================================

         call allocation (dt1, nppa(p), in_n(p), av_p(p), cl1(p), ca1(p)&
              &,cf1(p), sto_budg(:,p), day_storage(:,p), cl2(p), ca2(p)&
              &,cf2(p),litter_l(p),cwd(p),litter_fr(p),n_uptake(p), p_uptake(p)&
              &,lnr(:,p),end_pls)

         sto_budg(:,p) = day_storage(:,p)

         ! Se o PFT nao tem carbono goto 666-> TUDO ZERO
         if(end_pls) then
            no_pls_run = .true.
            dleaf(p) = 0.0
            dwood(p) = 0.0
            droot(p) = 0.0
            goto 666 ! gt hell
         endif

         if(ph(p) .eq. 0.0 .or. nppa(p) .eq. 0.0) then
            cue(p) = 0.0
         else
            cue(p) = nppa(p)/ph(p)
         endif

         dleaf(p) = cl2(p) - cl1(p)  !kg m-2
         dwood(p) = ca2(p) - ca1(p)
         droot(p) = cf2(p) - cf1(p)

!!!!====================================================
         !     Snow budget
         !     ===========
         smelt(p) = 2.63 + 2.55*temp + 0.0912*temp*prain !Snowmelt (mm/day)
         smelt(p) = amax1(smelt(p),0.)
         smelt(p) = amin1(smelt(p),s(p)+psnow)
         ds(p) = psnow - smelt(p)
         s(p) = s(p) + ds(p)

         !     Water budget
         !     ============
         if (ts.le.tice) then !Frozen soil
            g(p) = g(p) + w(p) !Soil moisture freezes
            w(p) = 0.0
            roff(p) = smelt(p) + prain !mm/day
            evap(p) = 0.0
            ph(p) = 0.0
            ar(p) = 0.0
            nppa(p) = 0.0
            laia(p) = 0.0
            rc2(p) = rcmin
            rm(p) = 0.0
            rg(p) = 0.0

         else                !Non-frozen soil
            w(p) = w(p) + g(p)
            g(p) = 0.0
            rimelt(p) = 0.0
            if (w(p).gt.wmax) then
               rimelt(p) = w(p) - wmax !Runoff due to soil ice melting
               w(p) = wmax
            endif

            roff(p) = runoff(w(p)/wmax)       !Soil moisture runoff (roff, mm/day)

            evap(p) = penman(p0,temp,rh,available_energy(temp),rc2(p)) !Actual evapotranspiration (evap, mm/day)
            dw(p) = prain + smelt(p) - evap(p) - roff(p)
            w(p) = w(p) + dw(p)
            if (w(p).gt.wmax) then
               roff(p) = roff(p) + (w(p) - wmax)
               w(p) = wmax
            endif
            if (w(p).lt.0.) w(p) = 0.
            roff(p) = roff(p) + rimelt(p) !Total runoff
         endif
!!!!====================================================
         litter_fr_aux = litter_fr(p) !* ocp_coeffs(p)
         litter_l_aux = litter_l(p) !* ocp_coeffs(p)
         cwd_aux = cwd(p) !* ocp_coeffs(p)

         ! if (real(n_uptake(p), r_4) .gt. in_n(p)) then
         !    print *, 'nupt gt inorg_n'
         !    call abort
         ! endif

         ! if (real(p_uptake(p), r_4) .gt. av_p(p)) then
         !    print *, 'pupt gt avail_p'
         !    call abort
         ! endif

         call carb3(ts, w(p)/wmax, litter_l_aux, cwd_aux, litter_fr_aux, real(lnr(:,p), r_4), clitter(:,p),&
                  & csoil(:, p), snr_internal(:,p),&
                  & av_p(p), in_n(p), in_p(p), so_p(p), litter_carbon_bdg(:,p),&
                  & soil_carbon_bdg(:,p), snr_aux(:,p), het_resp(p))



         !in_n(p) = in_n(p) - real(n_uptake(p), r_4)
         !av_p(p) = av_p(p) - real(p_uptake(p), r_4)

         do index = 1,8
            snr_internal(index,p) = snr_aux(index,p)
         enddo

         ! SIMULATE LOSS OF BIOMASS FROM POPULATION P DUE TO NEGATIVE NPP
         if(c_def(p) .gt. 0.0) then
            if(dt1(7) .gt. 0.0) then
               cl1_pft(p) = cl2(p) - ((c_def(p) * 1e-3) * 0.333333333)
               ca1_pft(p) = ca2(p) - ((c_def(p) * 1e-3) * 0.333333333)
               cf1_pft(p) = cf2(p) - ((c_def(p) * 1e-3) * 0.333333333)
            else
               cl1_pft(p) = cl2(p) - ((c_def(p) * 1e-3) * 0.5)
               ca1_pft(p) = 0.0
               cf1_pft(p) = cf2(p) - ((c_def(p) * 1e-3) * 0.5)
            endif
         else
            if(dt1(7) .gt. 0.0) then
               cl1_pft(p) = cl2(p)
               ca1_pft(p) = ca2(p)
               cf1_pft(p) = cf2(p)
            else
               cl1_pft(p) = cl2(p)
               ca1_pft(p) = 0.0
               cf1_pft(p) = cf2(p)
            endif
         endif

         ! FILL OUTPUT VARIABLES
         if (p .eq. 1) epavg = emax !mm/day
         w2(p)            = real(w(p),r_4)
         g2(p)            = real(g(p),r_4)
         s2(p)            = real(s(p),r_4)
         smavg(p)         = smelt(p)
         ruavg(p)         = roff(p)    ! mm day-1
         evavg(p)         = evap(p)    ! mm day-1
         phavg(p)         = ph(p)      !kgC/m2/day
         aravg(p)         = ar(p)      !kgC/m2/year
         nppavg(p)        = nppa(p)    !kgC/m2/day
         laiavg(p)        = laia(p)
         rcavg(p)         = rc2(p)     ! s m -1
         f5avg(p)         = f5(p)
         rmavg(p)         = rm(p)
         rgavg(p)         = rg(p)

         ! CVEG POOLS
         cleafavg_pft(p)  = cl1_pft(p)
         cawoodavg_pft(p) = ca1_pft(p)
         cfrootavg_pft(p) = cf1_pft(p)
         wueavg(p)        = wue(p)
         cueavg(p)        = cue(p)
         c_defavg(p)      = c_def(p)
         ocpavg(p)        = ocp_coeffs(p)

         ! CSOIL POOLS
         soilc(1:2,p)     = litter_carbon_bdg(:, p)
         soilc(3:4,p)     = soil_carbon_bdg(:, p)
         nupt(p)          = n_uptake(p)
         pupt(p)          = p_uptake(p)
         inorganic_p(p)   = in_p(p)
         inorganic_n(p)   = in_n(p)
         available_p(p)   = av_p(p)
         sorbed_p(p)      = so_p(p)

         do index = 1,8
            snr(index,p)  = snr_internal(index,p)
         enddo

         no_pls_run = .false.
666      continue
         if(no_pls_run) then
            ! All*** outputs are set to zero except epavg
            if(p .eq. 1) epavg = emax
            smavg(p) = 0.0
            ruavg(p) = 0.0
            evavg(p) = 0.0
            rcavg(p) = 0.0
            phavg(p) = 0.0
            aravg(p) = 0.0
            nppavg(p) = 0.0
            wueavg(p) = 0.0
            cueavg(p) = 0.0
            c_defavg(p) = 0.0
            laiavg(p) = 0.0
            f5avg(p) = 0.0
            rmavg(p) = 0.0
            rgavg(p) = 0.0
            vcmax(p) = 0.0
            specific_la(p) = 0.0
            sto_budg(:,p) =  0.0
            nupt(p) = 0.0
            pupt(p) = 0.0
            litter_l(p) = 0.0
            cwd(p) = 0.0
            litter_fr(p) = 0.0
            lnr(:,p) = 0.0
            ocpavg(p) = 0.0
            cleafavg_pft(p)  = 0.0
            cawoodavg_pft(p) = 0.0
            cfrootavg_pft(p) = 0.0
            cl1_pft(p) = 0.0
            ca1_pft(p) = 0.0
            cf1_pft(p) = 0.0
            w2(p) = 0.0
            g2(p) = 0.0
            s2(p) = 0.0
         endif

      enddo ! end pls_loop (p)

!      call pft_area_frac(cl1_pft, cf1_pft, ca1_pft, ocp_coeffs2, ocp_wood2)
      ! ocpavg = (ocp_coeffs + ocp_coeffs2

      ! CLEAN NANs OF some outputs
      do p = 1,npls
         if(isnan(cleafavg_pft(p))) cleafavg_pft(p) = 0.0
         if(isnan(cawoodavg_pft(p))) cawoodavg_pft(p) = 0.0
         if(isnan(cfrootavg_pft(p))) cfrootavg_pft(p) = 0.0
         if(cleafavg_pft(p) + 1 .eq. cleafavg_pft(p)) cleafavg_pft(p) = 0.0
         if(cawoodavg_pft(p) + 1 .eq. cawoodavg_pft(p)) cawoodavg_pft(p) = 0.0
         if(cfrootavg_pft(p) + 1 .eq. cfrootavg_pft(p)) cfrootavg_pft(p) = 0.0
      end do

      return
   end subroutine daily_budget

end module budget
