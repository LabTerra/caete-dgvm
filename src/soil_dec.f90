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
   use soil_pools
   implicit none

   ! Turnover Rates  == residence_time⁻¹ (years⁻¹)
   real(r_4), dimension(4) :: tr_c = (/10.0, 100.0, 1000.0, 1000.0/)
   !=========================================================================
   ! FUNCTIONS AND SUBROUTINES DEFINED IN SOIL_DEC MODULE
   public :: carbon3          ! Subroutine that calculates the C:N:P decay dynamics
   public :: carbon_decay     ! Carbon decay function in response to Temperarute
   public :: water_effect     ! Soil water content effect on C decay
   public :: sorbed_p_equil   ! Fucntion that caculates the equilibrium between Mineralized P and Sorbed P
   public :: bnf              ! Biological Nitrogen Fixation
   ! TODO
   ! public :: am_acivity       ! Arbuscular Micorrizal activity (Generate an output)
   ! public :: ptase_activity   ! Acivity of PTASE mediated by plant inputs of C
   ! public :: exudte_ activity ! Desorption of P from plant exudation of C
contains

   subroutine carbon3(tsoil, water_sat, leaf_litter, coarse_wd,&
                    &        root_litter, lnr, cl, cs, &
                    &        snr_in, avail_p, inorg_n, inorg_p,&
                    &        sorbed_p, avail_p_out, inorg_p_out,&
                    &        inorg_n_out, sorbed_p_out, cl_out, cs_out, snr, hr)

      real(r_4),parameter :: clit_atm = 0.7
      real(r_4),parameter :: cwd_atm = 0.22

      ! POOLS OF LITTER AND SOIL
      integer(i_4),parameter :: pl=2,ps=2
      integer(i_4) :: index, i

      !     Inputs
      !     ------
      real(r_4),intent(in) :: tsoil, water_sat       ! soil temperature (°C); soil water relative content (dimensionless)
      real(r_8),intent(in) :: leaf_litter
      real(r_8),intent(in) :: coarse_wd
      real(r_8),intent(in) :: root_litter            ! Mass of C comming from living pools g(C)m⁻²
      ! (lln2c),(rln2c),(cwdn2c),(llp2c),(rlp2c),(cwdp2c)
      real(r_4),dimension(6),intent(in) :: lnr       ! g(Nutrient) g(C)⁻¹ Incoming Nutrient Ratios

      real(r_4),dimension(pl),intent(in) :: cl       ! Litter carbon (gC/m2) State Variable -> The size of the carbon pools
      real(r_4),dimension(ps),intent(in) :: cs       ! Soil carbon (gC/m2)   State Variable -> The size of the carbon pools
     ! Soil Nutrient Ratios :: variables(8)            [(l1n2c),(l2n2c),(c1dn2c),(c2n2c),(l1p2c),(l2p2c),(c1p2c),(c2p2c)]
      real(r_4),dimension(8), intent(in) :: snr_in   ! Current soil nutrient ratios
      !real(r_4),intent(in) :: nupt, pupt             ! Nitrogen Uptake; Phosphorus Uptake (g m⁻² day⁻¹)

      !     Outputs
      !     -------
      real(r_4),intent(in) :: avail_p
      real(r_4),intent(in) :: inorg_p                  ! Pool of N biomineralized (gm⁻²)
      real(r_4),intent(in) :: inorg_n                  ! Pool of P biomineralized
      real(r_4),intent(in) :: sorbed_p                 ! Sorbed P - Secondary Mineral P

      real(r_4),intent(out) :: avail_p_out
      real(r_4),intent(out) :: inorg_p_out
      real(r_4),intent(out) :: inorg_n_out
      real(r_4),intent(out) :: sorbed_p_out

      real(r_4),dimension(pl),intent(out) :: cl_out       ! g(C)m⁻² State Variable -> The size of the carbon pools
      real(r_4),dimension(ps),intent(out) :: cs_out       ! State Variable -> The size of the carbon pools
      real(r_4),dimension(8), intent(out) :: snr          ! Updated Soil pools Nutrient to C ratios
      real(r_4),intent(out) :: hr                         ! Heterotrophic (microbial) respiration (gC/m2/day)

      real(r_4),dimension(4) :: nmass_org = 0.0
      real(r_4),dimension(4) :: pmass_org = 0.0

      !real(r_4),dimension(4) :: soil_nr    ! Soil pools Nutrient Ratios (2N + 2P)
      !real(r_4),dimension(4) :: litt_nr    ! litter pools Nutrient Ratios (2N + 2P)
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
      real(r_4) :: c_next_pool, n_next_pool, p_next_pool
      real(r_4) :: n_min_resp_lit, p_min_resp_lit
      real(r_4) :: incomming_c_lit, incomming_n_lit, incomming_p_lit
      real(r_4) :: update_c, update_n, update_p
      real(r_4) :: leaf_l, cwd, root_l, total_litter  ! Mass of C comming from living pools g(C)m⁻²

      ! START
      nutri_min_n = 0.0
      nutri_min_p = 0.0
      aux_ratio_n = 0.0
      aux_ratio_p = 0.0

      frac1 = 0.95
      frac2 = 1.0 - frac1


!     ! CARBON AND NUTRIENTS COMMING FROM VEGETATION

      leaf_l = real(leaf_litter, r_4)
      cwd    = real(coarse_wd, r_4)
      root_l = real(root_litter, r_4)

      ! Soil Nutrient Ratios structure:
      ! [ 1(l1n2c),2(l2n2c),3(c1dn2c),4(c2n2c),5(l1p2c),6(l2p2c),7(c1p2c),8(c2p2c)]

      ! find nutrient mass/area) : litter fluxes[ML⁻²] * litter nutrient ratios
      ! (lnr) [MM⁻¹]
      leaf_n  = leaf_l * lnr(1) ! g(nutrient) m-2
      froot_n = root_l * lnr(2)
      wood_n  = cwd    * lnr(3)
      leaf_p  = leaf_l * lnr(4)
      froot_p = root_l * lnr(5)
      wood_p  = cwd    * lnr(6)


      ! C:N:P CYCLING NUMERICAL SOLUTION
      ! ORGANIC NUTRIENTS in SOIL
      ! Soil Nutrient Ratios to output: Set to 0.0
      snr = 0.0

      ! Soil Nutrient ratios and organic nutrients g m-2
      do  i = 1,8
         snr_aux(i) = snr_in(i)
         if(isnan(snr_in(i))) snr_aux(i) = 0.0
         if(snr_in(i) .eq. snr_in(i) - 1.0) snr_aux(i) = 0.0
      enddo

      nmass_org(1) = snr_aux(1) * cl(1)                                      ! g(N)m-2
      nmass_org(2) = snr_aux(2) * cl(2)
      nmass_org(3) = snr_aux(3) * cs(1)
      nmass_org(4) = snr_aux(4) * cs(2)

      pmass_org(1) = snr_aux(5) * cl(1)                                      ! g(P)m-2
      pmass_org(2) = snr_aux(6) * cl(2)
      pmass_org(3) = snr_aux(7) * cs(1)
      pmass_org(4) = snr_aux(8) * cs(2)

      total_litter = leaf_l + root_l + cwd                                   ! g m-2

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

      !LITTER I
      cl_out(1) = cl(1) - cdec(1)                                          ! Update Litter Carbon 1

      ! Mineralization
      ! C
      ! Release od CO2
      het_resp(1) = cdec(1) * clit_atm                                       ! Heterotrophic respiration ! processed (dacayed) Carbon lost to ATM
      ! Carbon going to LITTER 2
      c_next_pool = cdec(1) - het_resp(1)                                    ! Carbon going to cl_out(2) (next pool)

      ! N
      ! N mineralized by the release of CO2
      n_min_resp_lit = het_resp(1) * snr_aux(1)
      ! N going to the LITTER II
      n_next_pool =  c_next_pool * snr_aux(1)

      ! UPDATE N in Organic MAtter Litter I
      nmass_org(1) = nmass_org(1) - (n_min_resp_lit + n_next_pool)

      ! UPDATE THE INORGANIC N POOL
      inorg_n_out = inorg_n + n_min_resp_lit
      n_min_resp_lit = 0.0

      ! P
      ! P mineralized by the release of CO2
      p_min_resp_lit = het_resp(1) * snr_aux(5)
      ! P going to the next pool
      p_next_pool = c_next_pool * snr_aux(5)

      ! UPDATE P in Organic MAtter LITTER I
      pmass_org(1) = pmass_org(1) - (p_min_resp_lit + p_next_pool)
      ! UPDATE THE INORGANIC N POOL
      inorg_p_out = inorg_p + p_min_resp_lit
      p_min_resp_lit = 0.0
      ! END OF MINERALIZATION PROCESS

      ! UPDATE CNP ORGANIC POOLS
      ! C
      incomming_c_lit = (frac1 * leaf_l) + (frac1 * root_l)                  ! INcoming Carbon from vegetation g m-2
      cl_out(1) = cl_out(1) + incomming_c_lit
      !N
      incomming_n_lit = (frac1 * leaf_n) + (frac1 * froot_n)                 ! Iincoming N
      nmass_org(1) = nmass_org(1) + incomming_n_lit
      !P
      incomming_p_lit = (frac1 * leaf_p) + (frac1 * froot_p)                 ! Incoming P
      pmass_org(1) = pmass_org(1) + incomming_p_lit

      ! UPDATE SNR
      snr(1) = nmass_org(1) / cl_out(1)
      snr(5) = pmass_org(1) / cl_out(1)
      ! END LITTER 1

      ! CLEAN AUX VARIABLES
      incomming_c_lit = 0.0
      incomming_n_lit = 0.0
      incomming_p_lit = 0.0
      update_c = c_next_pool
      c_next_pool = 0.0
      update_n = n_next_pool
      n_next_pool = 0.0
      update_p = p_next_pool
      p_next_pool = 0.0


      ! LITTER II
      cl_out(2) = cl(2) - cdec(2)

      ! Mineralization
      ! C
      ! Release od CO2
      het_resp(2) = cdec(2) * clit_atm                                       ! Heterotrophic respiration ! processed (dacayed) Carbon lost to ATM
      ! Carbon going to SOIL I
      c_next_pool = cdec(2) - het_resp(2)                                    ! Carbon going to cl_out(2) (next pool)

      ! N
      ! N mineralized by the release of CO2
      n_min_resp_lit = het_resp(2) * snr_aux(2)
      ! N going to the SOIL I
      n_next_pool =  c_next_pool * snr_aux(2)
      ! UPDATE N in Organic MAtter LITTER 2
      nmass_org(2) = nmass_org(2) - (n_min_resp_lit + n_next_pool)
      ! UPDATE THE INORGANIC N POOL
      inorg_n_out = inorg_n_out + n_min_resp_lit
      n_min_resp_lit = 0.0

      ! P
      ! P mineralized by the release of CO2
      p_min_resp_lit = het_resp(2) * snr_aux(6)
      ! P going to the next pool SOIL I
      p_next_pool = c_next_pool * snr_aux(6)

      ! UPDATE P ORGANIC POOL KITTER 2
      pmass_org(2) = pmass_org(2) - (p_min_resp_lit + p_next_pool)
      ! UPDATE THE INORGANIC N POOL
      inorg_p_out = inorg_p_out + p_min_resp_lit
      p_min_resp_lit = 0.0
      ! END OF MINERALIZATION PROCESS

            ! UPDATE CNP ORGANIC POOLS
      ! C
      incomming_c_lit = (frac2 * leaf_l) + (frac2 * root_l)                  ! INcoming Carbon from vegetation g m-2
      cl_out(2) = cl_out(2) + incomming_c_lit + update_c
      update_c = 0.0
      !N
      incomming_n_lit = (frac2 * leaf_n) + (frac2 * froot_n)                 ! Iincoming N
      nmass_org(2) = nmass_org(2) + incomming_n_lit + update_n
      update_n = 0.0
      !P
      incomming_p_lit = (frac2 * leaf_p) + (frac2 * froot_p)                 ! Incoming P
      pmass_org(2) = pmass_org(2) + incomming_p_lit + update_p
      update_p = 0.0

      ! UPDATE SNR
      snr(2) = nmass_org(2) / cl_out(2)
      snr(6) = pmass_org(2) / cl_out(2)
      ! END LITTER II

      ! CLEAN AUX VARIABLES
      incomming_c_lit = 0.0
      incomming_n_lit = 0.0
      incomming_p_lit = 0.0
      update_c = c_next_pool
      c_next_pool = 0.0
      update_n = n_next_pool
      n_next_pool = 0.0
      update_p = p_next_pool
      p_next_pool = 0.0

      !SOIL I The same steps commented for the litter pools

      cs_out(1) = cs(1) - cdec(3)

      ! Mineralization
      ! C
      ! Release od CO2
      het_resp(3) = cdec(3) * clit_atm                                       ! Heterotrophic respiration ! processed (dacayed) Carbon lost to ATM
      ! Carbon going to SOIL 2
      c_next_pool = cdec(3) - het_resp(3)

      ! N
      ! N mineralized by the release of CO2
      n_min_resp_lit = het_resp(3) * snr_aux(3)
      ! N going to the SOIL II
      n_next_pool =  c_next_pool * snr_aux(3)
      ! UPDATE N in Organic MAtter SOIL I
      nmass_org(3) = nmass_org(3) - (n_min_resp_lit + n_next_pool)
      ! UPDATE THE INORGANIC N POOL
      inorg_n_out = inorg_n_out + n_min_resp_lit
      n_min_resp_lit = 0.0

      ! P
      ! P mineralized by the release of CO2
      p_min_resp_lit = het_resp(3) * snr_aux(7)
      ! P going to the next pool
      p_next_pool = c_next_pool * snr_aux(7)
      pmass_org(3) = pmass_org(3) - (p_min_resp_lit + p_next_pool)
      ! UPDATE THE INORGANIC N POOL
      inorg_p_out = inorg_p_out + p_min_resp_lit
      p_min_resp_lit = 0.0
      ! END OF MINERALIZATION PROCESS

      ! UPDATE CNP ORGANIC POOLS
      ! C
      incomming_c_lit = 0.0                  ! INcoming Carbon from vegetation g m-2
      cs_out(1) = cs_out(1) + incomming_c_lit + update_c
      update_c = 0.0
      !N
      incomming_n_lit = 0.0                 ! Iincoming N
      nmass_org(3) = nmass_org(3) + incomming_n_lit + update_n
      update_n = 0.0
      !P
      incomming_p_lit = 0.0                 ! Incoming P
      pmass_org(3) = pmass_org(3) + incomming_p_lit + update_p
      update_p = 0.0

      ! UPDATE SNR
      snr(3) = nmass_org(3) / cs_out(1)
      snr(7) = pmass_org(3) / cs_out(1)
      ! END SOIL 1

      ! CLEAN AUX VARIABLES
      incomming_c_lit = 0.0
      incomming_n_lit = 0.0
      incomming_p_lit = 0.0
      update_c = c_next_pool
      c_next_pool = 0.0
      update_n = n_next_pool
      n_next_pool = 0.0
      update_p = p_next_pool
      p_next_pool = 0.0

      !SOIL II
      cs_out(2) = cs(2) - cdec(4)

      ! Mineralization
      ! C
      ! Release od CO2
      het_resp(4) = cdec(4)                                       ! Heterotrophic respiration ! processed (dacayed) Carbon lost to ATM
      ! Carbon going to SOIL 2

      ! N
      ! N mineralized by the release of CO2
      n_min_resp_lit = het_resp(4) * snr_aux(4)
      ! UPDATE N in Organic MAtter SOIL II
      nmass_org(4) = nmass_org(4) - n_min_resp_lit
      ! UPDATE THE INORGANIC N POOL
      inorg_n_out = inorg_n_out + n_min_resp_lit
      n_min_resp_lit = 0.0

      ! P
      ! P mineralized by the release of CO2
      p_min_resp_lit = het_resp(4) * snr_aux(8)
      ! UPDATE ORGANIC P POOL SOIL II
      pmass_org(4) = pmass_org(4) - p_min_resp_lit
      ! UPDATE THE INORGANIC N POOL
      inorg_p_out = inorg_p_out + p_min_resp_lit
      p_min_resp_lit = 0.0
      ! END OF MINERALIZATION PROCESS

      ! UPDATE CNP ORGANIC POOLS
      ! C
      incomming_c_lit = 0.0                  ! INcoming Carbon from vegetation g m-2
      cs_out(2) = cs_out(2) + incomming_c_lit + update_c
      update_c = 0.0
      !N
      incomming_n_lit = 0.0                 ! Iincoming N
      nmass_org(4) = nmass_org(4) + incomming_n_lit + update_n
      update_n = 0.0
      !P
      incomming_p_lit = 0.0                 ! Incoming P
      pmass_org(4) = pmass_org(4) + incomming_p_lit + update_p
      update_p = 0.0

      ! UPDATE SNR
      snr(4) = nmass_org(4) / cs_out(2)
      snr(8) = pmass_org(4) / cs_out(2)

      ! FINAL CALCULATIONS

      ! BNF
      inorg_n_out = inorg_n_out + bnf(0.0)

      ! INORGANIC P DYNAMICS

      sorbed_p_out = sorbed_p_equil(inorg_p_out)

      ! Include PTASE AND EXUDATES HERE
      avail_p_out = inorg_p_out - sorbed_p_out

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
