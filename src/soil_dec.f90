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

module soil_dec
   ! This module represents a global struct (This will only work if isolated in a sinlge process)
   ! To be used only in caete_dyn  and daily_budget subroutines
   use types
   use global_par

   implicit none
   private


   ! CONSTANTS FOR CARBON DECAY
   ! Turnover Rates  == residence_time⁻¹ (years⁻¹)
   real(r_4), dimension(4) :: tr_c = (/2.0, 35.0, 250.0, 5000.0/) !litter I (1) litter II (2) soilI (3) soil II (4)
!=========================================================================

   ! ## STATE VARIABLES
   ! TODO create getter and setter functions for these variables
   real(r_8),public,dimension(4) :: soc_soil_dec = 0.1D0
   real(r_8),public,dimension(4, npls) :: soil_carbon = 0.1D0 ! soil carbon pools
!=========================================================================
!
   ! These are global variables that are initialized in caete_init.
   ! New inputs -- Estimated from literature for Manaus Region
   ! For SPINUP only
   real(r_8),public :: available_n = 0.30775999D0 ! g m-2 Xu et al. 2013 ?
   real(r_8),public :: available_p = 0.204299955D0  ! g m-2 Yang et al., 2013

   real(r_8),public, dimension(npls) :: inorganic_n, inorganic_p, sorbed_p_pls
   real(r_8),public, dimension(npls) :: available_n_pls, available_p_pls
   real(r_8),public, dimension(npls) :: occp_coeffs

   real(r_8),public :: available_n_init = 0.30775999D0 ! g m-2 Xu et al. 2013 ?
   real(r_8),public :: available_p_init = 0.204299955D0  ! g m-2 Yang et al., 2013

   ! Internal Variables storing POOLS (IN)organic Nutrients
   real(r_8),public :: in_n
   real(r_8),public :: in_p
   real(r_8),public :: so_p
   real(r_8),public, dimension(npls) :: inorg_p = 0.0D0            ! Pool of N biomineralized (gm⁻²)
   real(r_8),public, dimension(npls) :: inorg_n = 0.0D0            ! Pool of P biomineralized
   real(r_8),public, dimension(npls) :: sorbed_p = 0.0D0           ! Sorbed P - Secondary Mineral P

   ! Available pools
!   real(r_4),private :: avail_p = 0.0            ! Avilable P
!   real(r_4),private :: avail_n = 0.0            ! Avilable N


   real(r_8),public,dimension(4,npls) :: nmass_org_pls = 0.0D0
   real(r_8),public,dimension(4,npls) :: pmass_org_pls = 0.0D0
   ! Soil nutrient Ratios
   real(r_8),public,dimension(4,npls) :: soil_nr    ! Soil pools Nutrient Ratios (2N + 2P)
   real(r_8),public,dimension(4,npls) :: litt_nr    ! litter pools Nutrient Ratios (2N + 2P)
   real(r_8),public,dimension(8,npls) :: snr_pls

   ! FUNCTIONS AND SUBROUTINES DEFINED IN SOIL_DEC MODULE
   public :: carbon3         ! Subroutine that calculates the C:N:P decay dynamics
   public :: carbon_decay    ! Carbon decay function in response to Temperarute
   public :: water_effect    ! Soil water content effect on C decay
   public :: sorbed_p_equil  ! Fucntion that caculates the equilibrium between Mineralized P and Sorbed P
   public :: bnf             ! Biological Nitrogen Fixation
   public :: cwm

contains


   subroutine carbon3(pls, ocp, tsoil, water_sat, leaf_l, cwd, root_l, lnr, nupt, pupt,&
                  & clitter, csoil, snr, hr)
      ! CARBON3 <- SOIL DECOMPOSITION MODEL FOR CAETÊ

      real(r_4),parameter :: clit_atm = 0.7
      real(r_4),parameter :: cwd_atm = 0.22

      ! POOLS OF LITTER AND SOIL
      integer(i_4),parameter :: pl=2,ps=2
      integer(i_4) :: index

      !     Variables
      !     =========
      !     Inputs
      !     ------
      integer(i_4),intent(in) ::           pls
      real(r_8),intent(in) ::              ocp
      real(r_4),intent(in) ::              tsoil
      real(r_4),intent(in) ::              water_sat       ! soil temperature (°C); soil water relative content (dimensionless)
      real(r_8),intent(in) ::              leaf_l
      real(r_8),intent(in) ::              cwd
      real(r_8),intent(in) ::              root_l    ! Mass of C comming from living pools g(C)m⁻²
      real(r_8),dimension(6),intent(in) :: lnr       ! g(Nutrient) g(C)⁻¹ Incoming Nutrient Ratios
      real(r_4),intent(in) ::              nupt
      real(r_4),intent(in) ::              pupt             ! Nitrogen Uptake; Phosphorus Uptake (g m⁻² day⁻¹)

      !     Outputs
      !     -------
      real(r_4),dimension(pl),intent(out) :: clitter       ! Litter carbon (gC/m2) State Variable -> The size of the carbon pools
      real(r_4),dimension(ps),intent(out) :: csoil       ! Soil carbon (gC/m2)   State Variable -> The size of the carbon pools
      real(r_8),dimension(8), intent(out) :: snr          ! Soil pools Nutrient to C ratios
      real(r_8),intent(out) ::               hr                    ! Heterotrophic (microbial) respiration (gC/m2/day)

      !TODO ! Insert output: Total mineralized N and P
      !INTERNAL VARIABLES
      real(r_8),dimension(pl) :: pl_nitrogen = 0.0   ! Nitrogênio da serapilheira -
      real(r_8),dimension(pl) :: pl_phosphorus = 0.0 ! FORFI do serapilheira
      real(r_8),dimension(ps) :: ps_nitrogen = 0.0   ! & so forth
      real(r_8),dimension(ps) :: ps_phosphorus = 0.0

      real(r_8) :: leaf_n  ! Mass of Nutrients in MO comming from Vegetation  g(Nutrient) m⁻²
      real(r_8) :: froot_n
      real(r_8) :: wood_n
      real(r_8) :: leaf_p
      real(r_8) :: froot_p
      real(r_8) :: wood_p

      real(r_8) :: water_modifier ! Multiplicator for water influence on C decay

      real(r_8) :: frac1,frac2    ! Constants for litter partitioning

      real(r_8),dimension(pl+ps) :: het_resp, cdec
      real(r_8),dimension(pl+ps) :: aux_ratio_n, aux_ratio_p
      real(r_8),dimension(pl+ps) :: nutri_min_n, nutri_min_p

      real(r_8),dimension(pl) :: cl = 0.0D0
      real(r_8),dimension(ps) :: cs = 0.0D0

      !Auxiliary variables
      real(r_8) :: aux1, aux2, aux3, aux4, auxn, auxp
      real(r_8),dimension(4) :: nmass_org, pmass_org
      ! START

      occp_coeffs(pls) = ocp
      ! Must update these variables at the end of execution
      cl = soil_carbon(1:2, pls)
      cs = soil_carbon(3:4, pls)
      nmass_org = nmass_org_pls(:, pls)
      pmass_org = pmass_org_pls(:, pls)


      ! Start some variables
      nutri_min_n = 0.0
      nutri_min_p = 0.0
      aux_ratio_n = 0.0
      aux_ratio_p = 0.0

      frac1 = 0.8
      frac2 = 1.0 - frac1

      ! Soil Nutrient Ratios: Set to 0.0
      snr = 0.0

      ! find nutrient mass/area) : litter fluxes[ML⁻²] * litter nutrient ratios
      ! (lnr) [MM⁻¹]
      leaf_n    = real(leaf_l * lnr(1), r_8) ! g(nutrient) m-2
      froot_n   = real(root_l * lnr(2), r_8)
      wood_n    = real(cwd * lnr(3), r_8)
      leaf_p    = real(leaf_l * lnr(4), r_8)
      froot_p   = real(root_l * lnr(5), r_8)
      wood_p    = real(cwd * lnr(6), r_8)

      ! FIRST OF ALL calculate dacay from pools
      ! CARBON DECAY
      water_modifier = water_effect(water_sat)
      do index = 1,4
         if(index .lt. 3) then
            ! FOR THE 2 LITTER POOLS
            cdec(index) = carbon_decay(q10,tsoil,cl(index),tr_c(index)) * water_modifier
         else
            ! FOR THE 3 CARBON POOLS
            cdec(index) = carbon_decay(q10,tsoil,cs(index-2),tr_c(index)) * water_modifier
         endif
      enddo

      ! C:N:P CYCLING

      !LITTER I
      aux1 = (frac1 * leaf_l) + (frac1 * root_l)                             ! INcoming Carbon from vegetation
      aux2 = cdec(1) * clit_atm                                              ! processed (dacayed) Carbon lost to ATM
      aux3 = cdec(1) - aux2                                                  ! Carbon going to cl_out(2) (next pool)
      cl(1) = cl(1) - real(cdec(1) + aux1, r_4)                                   ! Update Litter Carbon 1
      het_resp(1) = real(aux2, r_4)                                                     ! Heterotrophic respiration
      aux4 = aux3                                                            ! Carbon going to the next pool
      pl_nitrogen(1) = nmass_org(1) + (leaf_n * frac1) + (froot_n * frac1)   ! Calculate the organic N/P pool
      pl_phosphorus(1) = pmass_org(1) + (leaf_p * frac1) + (froot_p * frac1) ! g(N/P)m-2
      nmass_org(1) = pl_nitrogen(1)                                          ! Update Organic N pool
      pmass_org(1) = pl_phosphorus(1)                                        ! Update Organic P pool

      ! LITTER II
      aux1 = (frac2 * leaf_l) + (frac2 * root_l) + (cwd * frac2) + aux4         ! Incoming Carbon
      aux2 = cdec(2) * clit_atm                                                 ! C to ATM
      aux3 = cdec(2) - aux2                                                     ! C going to csoil(1)
      cl(2) = cl(2) - real(cdec(2) + aux1, r_4)                                      ! Update Litter Carbon 2
      het_resp(2) = real(aux2, r_4)                                                        ! Heterotrophic respiration
      aux4 = aux3                                                               ! Carbon going to the next pool
      pl_nitrogen(2) = nmass_org(2) + (leaf_n * frac2) + (froot_n * frac2)&     ! Calculate the organic N/P pool
                        & + (wood_n *frac2)
      pl_phosphorus(2) = pmass_org(2) + (leaf_p * frac2) + (froot_p * frac2)&   ! g(N/P)m-2
                        & + (wood_p*frac2)
      nmass_org(2) = pl_nitrogen(2)                                             ! Update Organic N pool
      pmass_org(2) = pl_phosphorus(2)                                           ! Update Organic P pool

      !SOIL I The same steps commented for the litter pools
      aux1 = (frac1 * cwd) + aux4                             ! Incoming Carbon
      aux2 = cdec(3) * cwd_atm                                ! C to ATM
      aux3 = cdec(3) - aux2                                   ! C to csoil(2)
      cs(1) = cs(1) - real(cdec(3) + aux1, r_4)                    ! Update Soil Carbon 1
      het_resp(3) = real(aux2, r_4)                                      ! Het resp
      aux4 = aux3                                             ! Carbon going to the next pool
      ps_nitrogen(1) = nmass_org(3) + (wood_n * frac1)        ! Calculate the organic N/P pool
      ps_phosphorus(1) = pmass_org(3) + (wood_p * frac1)      ! g(N/P)m-2

      !Before update Organic pools - Transfer some nutrients to the last pool (2%)
      aux1 = 0.02D0 * ps_nitrogen(1)
      nmass_org(3) = ps_nitrogen(1) - aux1                     ! Update Organic N pool
      aux2 = 0.02D0 * ps_phosphorus(1)
      pmass_org(3) = ps_phosphorus(1) - aux2                         ! Update Organic P pool


      !SOIL II
      het_resp(4) = cdec(4) * cwd_atm            ! C to atm
      cs(2) = cs(2) - real(het_resp(4) + aux4, r_4)     ! Aux 4 is the carbon comming from preceding pool

      ps_nitrogen(2) = nmass_org(4) + aux1
      ps_phosphorus(2) = pmass_org(4) + aux2

      nmass_org(4) = ps_nitrogen(2)
      pmass_org(4) = ps_phosphorus(2)

      ! THIS SECTION - Mineralized nutrients =====================================================
      ! The amounts of Minerilized nutrients are dependent on N:C and P:C mass ratios of soil pools
      ! NUTRIENT RATIOS in SOIL
      do index = 1,4
         if(index .lt. 3) then
            aux_ratio_n(index) = nmass_org(index) / real(cl(index), r_8) ! g(N)g(C)-1
            aux_ratio_p(index) = pmass_org(index) / real(cl(index), r_8) ! g(P)g(C)-1
         else
            if (cs(index-2) .le. 0.0) then
               aux_ratio_n(index) = 0.0
            else
               aux_ratio_n(index) = nmass_org(index) / cs(index-2) ! g(N)g(C)-1
            endif
            if (cs(index-2) .le. 0.0) then
               aux_ratio_p(index) = 0.0
            else
                aux_ratio_p(index) = pmass_org(index) / cs(index-2) ! g(P)g(C)-1
            endif
         endif
      enddo

      ! **** To GLOBAL VARS
      litt_nr(1, pls) = aux_ratio_n(1)
      litt_nr(2, pls) = aux_ratio_n(2)
      soil_nr(1, pls) = aux_ratio_n(3)
      soil_nr(2, pls) = aux_ratio_n(4)
      litt_nr(3, pls) = aux_ratio_p(1)
      litt_nr(4, pls) = aux_ratio_p(2)
      soil_nr(3, pls) = aux_ratio_p(3)
      soil_nr(4, pls) = aux_ratio_p(4)

      ! OUTPUT SOIL NUTRIENT RATIOS
      do index = 1,8
         if (index .lt. 5) snr_pls(index,pls) = aux_ratio_n(index)
         if (index .gt. 4) snr_pls(index,pls) = aux_ratio_p(index-4)
      end do

      ! USE NUTRIENT RATIOS AND HET_RESP TO CALCULATE MINERALIZED NUTRIENTS
      do index=1,4
         nutri_min_n(index) = het_resp(index)*aux_ratio_n(index)
         nutri_min_p(index) = het_resp(index)*aux_ratio_p(index)
      enddo

      ! UPDATE N and P ORGANIC in SOIL POOLS
      do index=1,4
         nmass_org(index) = nmass_org(index) - nutri_min_n(index)
         pmass_org(index) = pmass_org(index) - nutri_min_p(index)
      enddo

      soil_carbon(1:2, pls) = cl
      soil_carbon(3:4, pls) = cs
      nmass_org_pls(:, pls) = nmass_org
      pmass_org_pls(:, pls) = pmass_org

      clitter(1) = real(cl(1), r_4)
      clitter(2) = real(cl(2), r_4)
      csoil(1) = real(cs(1), r_4)
      csoil(2) = real(cs(2), r_4)

      ! UPDATE INORGANIC POOLS
      inorg_n(pls) =  inorg_n(pls) + sum(nutri_min_n)
      ! Update available N pool
      ! BNF
      auxn = inorg_n(pls) - (nupt * 1.0) !+ auxn ! Global Variable
      ! INLCUDE SORPTION DYNAMICS
      inorg_p(pls) = inorg_p(pls) + sum(nutri_min_p)
      sorbed_p(pls) = sorbed_p_equil(inorg_p(pls))
      auxp = inorg_p(pls) - sorbed_p(pls) - (pupt * 1.0) !+ auxp

      if (inorg_p(pls) .lt. 0.0) inorg_p(pls) = 0.0
      if (inorg_n(pls) .lt. 0.0) inorg_n(pls) = 0.0

      if (auxp .lt. 0.0) auxp = 0.0
      if (auxn .lt. 0.0) auxn = 0.0

      available_n_pls(pls)    = auxn
      available_p_pls(pls)    = auxp
      inorganic_n(pls)        = inorg_n(pls)
      inorganic_p(pls)        = inorg_p(pls)
      sorbed_p_pls(pls)       = sorbed_p(pls)

      if(pls .eq. npls) then
         in_n         = cwm(inorganic_n, occp_coeffs)
         in_p         = cwm(inorganic_p, occp_coeffs)
         so_p         = cwm(sorbed_p_pls, occp_coeffs)
         available_n  = cwm(available_n_pls, occp_coeffs)
         available_p  = cwm(available_p_pls, occp_coeffs)
         snr          = snr_pls(:, pls)
         inorg_n(:)   = in_n
         inorg_p(:)   = in_p
         sorbed_p(:)  = so_p
         soc_soil_dec(1)     = cwm(real(soil_carbon(1,:),r_8), occp_coeffs)
         soc_soil_dec(2)     = cwm(real(soil_carbon(2,:),r_8), occp_coeffs)
         soc_soil_dec(3)     = cwm(real(soil_carbon(3,:),r_8), occp_coeffs)
         soc_soil_dec(4)     = cwm(real(soil_carbon(4,:),r_8), occp_coeffs)
      endif

      !OUTPUT  Heterotrophic respiration
      hr = sum(het_resp)

      ! FEATURES TO BE IMPLEMENTED

      ! P biochemical mineralization +
      ! P release "de-sorption" +

   end subroutine carbon3

   ! -------------------------------------------------------------

   function carbon_decay(q10_in,tsoil,c,residence_time) result(decay)
      !Based on carbon decay implemented in JeDi and JSBACH - Pavlick et al. 2012
         real(r_4),intent(in) :: q10_in           ! constant ~1.4
         real(r_4),intent(in) :: tsoil            ! Soil temperature °C
         real(r_8),intent(in) :: c                ! Carbon content per area g(C)m-2
         real(r_4),intent(in) :: residence_time   ! Pool turnover rate
         real(r_4) :: decay                       ! ML⁻²

         if(c .le. 0.0) then
            decay = 0.0
            return
         endif

         decay = (q10_in ** ((tsoil - 20.0) / 10.0)) * (real(c, r_4) / residence_time)
      end function carbon_decay


      function water_effect(theta) result(retval)
         ! Implement the Moyano function based on soil water content. Moyano et al. 2012;2013
         ! Based on the implementation of Sierra et al. 2012 (SoilR)
         ! This fucntion is ideal and was parametrized for low carbon soils

         real(r_4),intent(in) :: theta  ! Volumetric soil water content (cm³ cm⁻³)
         real(r_4),parameter :: k_a = 3.11, k_b = 2.42
         real(r_4) :: retval

         retval= (k_a * theta) - (k_b * theta**2)

      end function water_effect


      function sorbed_p_equil(arg) result(retval)
         ! Linear equilibrium between inorganic P and available P pool

         real(r_8), intent(in) :: arg
         real(r_8) :: retval

         retval = arg * ks
      end function sorbed_p_equil


      function bnf(c_amount) result(n_amount)

         real(r_4), intent(in) :: c_amount
         real(r_4) :: n_amount

         n_amount = c_amount * (1.0/30.0)

      end function bnf


      function cwm(var_arr, area_arr) result(retval)

         use types
         use global_par

         real(kind=r_8), dimension(npls), intent(in) :: var_arr, area_arr
         real(kind=r_8) :: retval
         retval = sum(var_arr * area_arr, mask = .not. isnan(var_arr))

      end function cwm

end module soil_dec
