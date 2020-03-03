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

   use types
   use global_par

   implicit none
   private


   ! CONSTANTS FOR CARBON DECAY
   ! Turnover Rates  == residence_time⁻¹ (years⁻¹)
   real(r_4), dimension(4) :: tr_c = (/2.0, 15.0, 150.0, 3000.0/)
   !litter I (1) litter II (2) soilI (3) soil II (4)
   real(r_4), public :: sp_available_p = 0.204299955    ! g m-2 Yang et al., 2013
   real(r_4), public :: sp_available_n = 0.30775999     ! g m-2 Xu et al. 2013 ?
   !=========================================================================
   ! FUNCTIONS AND SUBROUTINES DEFINED IN SOIL_DEC MODULE
   public :: carbon3         ! Subroutine that calculates the C:N:P decay dynamics
   public :: carbon_decay    ! Carbon decay function in response to Temperarute
   public :: water_effect    ! Soil water content effect on C decay
   public :: sorbed_p_equil  ! Fucntion that caculates the equilibrium between Mineralized P and Sorbed P
   public :: bnf             ! Biological Nitrogen Fixation

contains

   subroutine carbon3(tsoil, water_sat, leaf_litter, coarse_wd, root_litter, lnr, cl, cs, &
                    &  snr_in, nupt, pupt, avail_p, inorg_n, inorg_p,&
                    & sorbed_p, cl_out, cs_out, snr, hr)

      real(r_4),parameter :: clit_atm = 0.7
      real(r_4),parameter :: cwd_atm = 0.22

      ! POOLS OF LITTER AND SOIL
      integer(i_4),parameter :: pl=2,ps=2
      integer(i_4) :: index, i, j

      !     Inputs
      !     ------
      real(r_4),intent(in) :: tsoil, water_sat       ! soil temperature (°C); soil water relative content (dimensionless)
      real(r_8),intent(in) :: leaf_litter
      real(r_8),intent(in) :: coarse_wd
      real(r_8),intent(in) :: root_litter            ! Mass of C comming from living pools g(C)m⁻²
      real(r_4),dimension(6),intent(in) :: lnr       ! g(Nutrient) g(C)⁻¹ Incoming Nutrient Ratios

      real(r_4),dimension(pl),intent(in) :: cl       ! Litter carbon (gC/m2) State Variable -> The size of the carbon pools
      real(r_4),dimension(ps),intent(in) :: cs       ! Soil carbon (gC/m2)   State Variable -> The size of the carbon pools
      real(r_4),dimension(8), intent(in) :: snr_in   ! Current soil nutrient ratios
      real(r_4),intent(in) :: nupt, pupt             ! Nitrogen Uptake; Phosphorus Uptake (g m⁻² day⁻¹)

      !     Outputs
      !     -------
      real(r_4),intent(inout) :: avail_p
      real(r_4),intent(inout) :: inorg_p                  ! Pool of N biomineralized (gm⁻²)
      real(r_4),intent(inout) :: inorg_n                  ! Pool of P biomineralized
      real(r_4),intent(inout) :: sorbed_p                 ! Sorbed P - Secondary Mineral P
      real(r_4),dimension(pl),intent(out) :: cl_out       ! g(C)m⁻² State Variable -> The size of the carbon pools
      real(r_4),dimension(ps),intent(out) :: cs_out       ! State Variable -> The size of the carbon pools
      real(r_4),dimension(8), intent(out) :: snr          ! Updated Soil pools Nutrient to C ratios
      real(r_4),intent(out) :: hr                         ! Heterotrophic (microbial) respiration (gC/m2/day)

      real(r_4),dimension(4) :: nmass_org = 0.0
      real(r_4),dimension(4) :: pmass_org = 0.0

      real(r_4),dimension(4) :: soil_nr    ! Soil pools Nutrient Ratios (2N + 2P)
      real(r_4),dimension(4) :: litt_nr    ! litter pools Nutrient Ratios (2N + 2P)
      !TODO ! Insert output: Total mineralized N and P
      !INTERNAL VARIABLES
      real(r_4),dimension(pl) :: pl_nitrogen = 0.0   ! Nitrogênio da serapilheira -
      real(r_4),dimension(pl) :: pl_phosphorus = 0.0 ! FORFI do serapilheira
      real(r_4),dimension(ps) :: ps_nitrogen = 0.0   ! & so forth
      real(r_4),dimension(ps) :: ps_phosphorus = 0.0

      real(r_4) :: leaf_n  ! Mass of Nutrients in MO comming from Vegetation  g(Nutrient) m⁻²
      real(r_4) :: froot_n
      real(r_4) :: wood_n
      real(r_4) :: leaf_p
      real(r_4) :: froot_p
      real(r_4) :: wood_p

      real(r_4) :: water_modifier ! Multiplicator for water influence on C decay

      real(r_4) :: frac1,frac2    ! Constants for litter partitioning

      real(r_4),dimension(pl+ps) :: het_resp, cdec
      real(r_4),dimension(pl+ps) :: aux_ratio_n, aux_ratio_p
      real(r_4),dimension(pl+ps) :: nutri_min_n, nutri_min_p
      real(r_4),dimension(8) :: snr_aux

      !Auxiliary variables
      real(r_4) :: aux1, aux2, aux3, aux4
      real(r_4) :: leaf_l, cwd, root_l    ! Mass of C comming from living pools g(C)m⁻²


      ! START
      leaf_l = real(leaf_litter, r_4)
      cwd    = real(coarse_wd, r_4)
      root_l = real(root_litter, r_4)

      nutri_min_n = 0.0
      nutri_min_p = 0.0
      aux_ratio_n = 0.0
      aux_ratio_p = 0.0

      frac1 = 0.95
      frac2 = 1.0 - frac1

      ! Soil Nutrient Ratios structure:
      ! [ 1(l1n2c),2(l2n2c),3(c1dn2c),4(c2n2c),5(l1p2c),6(l2p2c),7(c1p2c),8(c2p2c)]

      ! Soil Nutrient Ratios to output: Set to 0.0
      snr = 0.0

      ! Soil Nutrient ratios and organic nutrients g m-2
      do  i = 1,8
         snr_aux(i) = snr_in(i)
         if(isnan(snr_in(i))) snr_aux(i) = 0.0D0
         if(snr_in(i) .eq. snr_in(i) - 1.0D0) snr_aux(i) = 0.0D0
      enddo

      nmass_org(1) = snr_aux(1) * cl(1)
      nmass_org(2) = snr_aux(2) * cl(2)
      nmass_org(3) = snr_aux(3) * cs(1)
      nmass_org(4) = snr_aux(4) * cs(2)

      pmass_org(1) = snr_aux(5) * cl(1)
      pmass_org(2) = snr_aux(6) * cl(2)
      pmass_org(3) = snr_aux(7) * cs(1)
      pmass_org(4) = snr_aux(8) * cs(2)


      ! find nutrient mass/area) : litter fluxes[ML⁻²] * litter nutrient ratios
      ! (lnr) [MM⁻¹]
      leaf_n  = leaf_l * lnr(1) ! g(nutrient) m-2
      froot_n = root_l * lnr(2)
      wood_n  = cwd    * lnr(3)
      leaf_p  = leaf_l * lnr(4)
      froot_p = root_l * lnr(5)
      wood_p  = cwd    * lnr(6)

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
      cl_out(1) = (cl(1) - cdec(1)) + aux1                                   ! Update Litter Carbon 1
      het_resp(1) = aux2                                                     ! Heterotrophic respiration
      aux4 = aux3                                                            ! Carbon going to the next pool
      pl_nitrogen(1) = nmass_org(1) + (leaf_n * frac1) + (froot_n * frac1)   ! Calculate the organic N/P pool
      pl_phosphorus(1) = pmass_org(1) + (leaf_p * frac1) + (froot_p * frac1) ! g(N/P)m-2
      nmass_org(1) = pl_nitrogen(1)                                          ! Update Organic N pool
      pmass_org(1) = pl_phosphorus(1)                                        ! Update Organic P pool

      ! LITTER II
      aux1 = (frac2 * leaf_l) + (frac2 * root_l) + (cwd * frac2) + aux4         ! Incoming Carbon
      aux2 = cdec(2) * clit_atm                                                 ! C to ATM
      aux3 = cdec(2) - aux2                                                     ! C going to csoil(1)
      cl_out(2) = (cl(2) - cdec(2)) + aux1                                      ! Update Litter Carbon 2
      het_resp(2) = aux2                                                        ! Heterotrophic respiration
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
      cs_out(1) = (cs(1) - cdec(3)) + aux1                    ! Update Soil Carbon 1
      het_resp(3) = aux2                                      ! Het resp
      aux4 = aux3                                             ! Carbon going to the next pool
      ps_nitrogen(1) = nmass_org(3) + (wood_n * frac1)        ! Calculate the organic N/P pool
      ps_phosphorus(1) = pmass_org(3) + (wood_p * frac1)      ! g(N/P)m-2

      !Before update Organic pools - Transfer some nutrients to the last pool (2%)
      aux1 = 0.02 * ps_nitrogen(1)
      nmass_org(3) = ps_nitrogen(1) - aux1                     ! Update Organic N pool
      aux2 = 0.02 * ps_phosphorus(1)
      pmass_org(3) = ps_phosphorus(1) - aux2                         ! Update Organic P pool


      !SOIL II
      het_resp(4) = cdec(4) * cwd_atm            ! C to atm
      cs_out(2) = cs(2) - het_resp(4) + aux4     ! Aux 4 is the carbon comming from preceding pool

      ps_nitrogen(2) = nmass_org(4) + aux1
      ps_phosphorus(2) = pmass_org(4) + aux2

      nmass_org(4) = ps_nitrogen(2)
      pmass_org(4) = ps_phosphorus(2)

      ! THIS SECTION - Mineralized nutrients =====================================================
      ! The amounts of Minerilized nutrients are dependent on N:C and P:C mass ratios of soil pools
      ! NUTRIENT RATIOS in SOIL
      do index = 1,4
         if(index .lt. 3) then
            aux_ratio_n(index) = nmass_org(index) / cl_out(index) ! g(N)g(C)-1
            aux_ratio_p(index) = pmass_org(index) / cl_out(index) ! g(P)g(C)-1
         else
            if (cs_out(index-2) .le. 0.0) then
               aux_ratio_n(index) = 0.0
            else
               aux_ratio_n(index) = nmass_org(index) / cs_out(index-2) ! g(N)g(C)-1
            endif
            if (cs_out(index-2) .le. 0.0) then
               aux_ratio_p(index) = 0.0
            else
                aux_ratio_p(index) = pmass_org(index) / cs_out(index-2) ! g(P)g(C)-1
            endif
         endif
      enddo

      ! **** To GLOBAL VARS
      litt_nr(1) = aux_ratio_n(1)
      litt_nr(2) = aux_ratio_n(2)
      soil_nr(1) = aux_ratio_n(3)
      soil_nr(2) = aux_ratio_n(4)
      litt_nr(3) = aux_ratio_p(1)
      litt_nr(4) = aux_ratio_p(2)
      soil_nr(3) = aux_ratio_p(3)
      soil_nr(4) = aux_ratio_p(4)

      ! OUTPUT SOIL NUTRIENT RATIOS
      do index = 1,8
         if (index .lt. 5) snr(index) = aux_ratio_n(index)
         if (index .gt. 4) snr(index) = aux_ratio_p(index-4)
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

      ! UPDATE INORGANIC POOLS
      inorg_n =  inorg_n + sum(nutri_min_n) - nupt
      ! BNF
      ! INLCUDE SORPTION DYNAMICS
      inorg_p = inorg_p + sum(nutri_min_p)
      sorbed_p = sorbed_p_equil(inorg_p)
      avail_p = inorg_p - sorbed_p - pupt !+ available_p

      if (inorg_p .lt. 0.0) inorg_p = 0.0
      if (inorg_n .lt. 0.0) inorg_n = 0.0
      if (avail_p .lt. 0.0) avail_p = 0.0

      hr = sum(het_resp)

   end subroutine carbon3


   function carbon_decay(q10_in,tsoil,c,residence_time) result(decay)
   !Based on carbon decay implemented in JeDi and JSBACH - Pavlick et al. 2012
      real(r_4),intent(in) :: q10_in           ! constant ~1.4
      real(r_4),intent(in) :: tsoil            ! Soil temperature °C
      real(r_4),intent(in) :: c                ! Carbon content per area g(C)m-2
      real(r_4),intent(in) :: residence_time   ! Pool turnover rate
      real(r_4) :: decay                       ! ML⁻²

      if(c .le. 0.0) then
         decay = 0.0
         return
      endif

      decay = (q10_in ** ((tsoil - 20.0) / 10.0)) * (c / residence_time)
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

      real(r_4), intent(in) :: arg
      real(r_4) :: retval

      retval = arg * ks
   end function sorbed_p_equil


   function bnf(c_amount) result(n_amount)

      real(r_4), intent(in) :: c_amount
      real(r_4) :: n_amount

      n_amount = c_amount * (1.0/30.0)

   end function bnf

end module soil_dec
