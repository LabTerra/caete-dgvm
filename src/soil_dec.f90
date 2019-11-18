! Author: João Paulo Darela FIlho- joao.darela@unesp.br


module soil_dec
   ! This module represents a global struct (This will only work if isolated in a sinlge process)
   ! To be used only in caete_dyn  and daily_budget subroutines
   use types
   use global_par

   implicit none
   private

   ! ## STATE VARIABLES TODO
   ! TODO create getter and setter functions for these variables
   real(r_4),public,dimension(2) :: litter_carbon = 0.0 ! Litter carbon pools
   real(r_4),public,dimension(2) :: soil_carbon = 0.0   ! soil carbon pools

!
   ! These are global variables that are initialized in caete_init.
   ! New inputs -- Estimated from literature for Manaus Region
   real(r_8),public :: available_n = 3.775999d-4 ! kg m-2 Xu et al. 2013 ?
   real(r_8),public :: available_p = 2.4299955d-4  ! kg m-2 Yang et al., 2013


   ! Internal Variables storing POOLS INorganic Nutrients
   real(r_4),private :: inorg_p
   real(r_4),private :: inorg_n

   real(r_4),private,dimension(4) :: nmass_org = 0.0
   real(r_4),private,dimension(4) :: pmass_org = 0.0
   ! Soil nutrient Ratios
   real(r_4),private,dimension(4) :: soil_nr    ! Soil pools Nutrient Ratios (2N + 2P)
   real(r_4),private,dimension(4) :: litt_nr    ! litter pools Nutrient Ratios (2N + 2P)

   !Auxiliary variables
   real(r_4),private :: aux1, aux2, aux3, aux4


! FUNCTIONS AND SUBROUTINES DEFINED IN SOIL_DEC MODULE
   public :: carbon3
   public :: scarbon_decayment
   public :: water_effect

   public :: get_inorgp, set_inorgp
   public :: get_inorgn, set_inorgn


contains
   ! GETTERS AND SETTERS FOR P AND N INORGANIC POOLS
   function get_inorgp() result(retval)

      real(r_4) :: retval

      retval = inorg_p

   end function get_inorgp

   function get_inorgn() result(retval)

      real(r_4) :: retval

      retval = inorg_n

   end function get_inorgn

   subroutine set_inorgp(arg)

      real(r_4),intent(in) :: arg

      inorg_p = arg

   end subroutine set_inorgp

   subroutine set_inorgn(arg)

      real(r_4),intent(in) :: arg

      inorg_n = arg

   end subroutine set_inorgn
 ! GETTERS AND SETTERS FOR P AND N INORGANIC POOLS
! -------------------------------------------------------------

   function scarbon_decayment(q10_in,tsoil,c,residence_time) result(decay)
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

      decay = (q10_in**((tsoil - 20.0) / 10.0)) * (c/residence_time)
   end function scarbon_decayment

   function water_effect(theta) result(retval)
   ! Implement the Moyano function based on soil water content. Moyano et al. 2013
   ! Based on the implementation of Sierra et al. 2012 (SoilR)
      real(r_4),intent(in) :: theta  ! Volumetric soil water content (cm³ cm⁻³)
      real(r_4),parameter :: k_a = 3.11, k_b = 2.42
      real(r_4) :: retval

      retval= (k_a * theta) - (k_b * theta**2)
   end function water_effect


   subroutine carbon3(tsoil,water_sat,leaf_l, cwd, root_l, lnr, cl, cs, &
                    &  cl_out, cs_out, snr, hr)
      ! From Pavlick et al. 2012 - JeDi/ Raddatz et al. 2007 - JSBACH

       ! Fraction of C that is broken by soil microbial activity and goes
      ! to Atmosphere (A.K.A. Heterothrophic respiration )
      real(r_4),parameter :: clit_atm = 0.80
      real(r_4),parameter :: cwd_atm = 0.40

      integer(i_4),parameter :: pl=2,ps=2
      integer(i_4) :: index

      !     Variables
      !     =========
      !     Inputs
      !     ------
      real(r_4),intent(in) :: tsoil, water_sat       ! soil temperature (°C); soil water relative content (dimensionless)
      real(r_4),intent(in) :: leaf_l                 ! Mass of C comming from living pools g(C)m⁻²
      real(r_4),intent(in) :: cwd                    ! Mass of C comming from living pools g(C)m⁻²
      real(r_4),intent(in) :: root_l                 ! Mass of C comming from living pools g(C)m⁻²
      real(r_4),dimension(6),intent(in) :: lnr       !g(Nutrient) g(C)⁻¹ Incoming Nutrient Ratios

      real(r_4),dimension(pl),intent(in) :: cl       !Litter carbon (gC/m2) State Variable -> The size of the carbon pools
      real(r_4),dimension(ps),intent(in) :: cs       !Soil carbon (gC/m2)   State Variable -> The size of the carbon pools

      !     Outputs
      !     -------
      real(r_4),dimension(pl),intent(out) :: cl_out  ! g(C)m⁻² State Variable -> The size of the carbon pools
      real(r_4),dimension(ps),intent(out) :: cs_out  !         State Variable -> The size of the carbon pools
      real(r_4),dimension(8), intent(out) :: snr ! Soil pools Nutrient to C ratios
      real(r_4),intent(out) :: hr                    !Heterotrophic (microbial) respiration (gC/m2/day)


      !TODO ! Insert output: Total mineralized N and P
      !INTERNAL VARIABLES
      real(r_4),dimension(pl+ps) :: tr_c
      real(r_4),dimension(pl) :: pl_nitrogen = 0.0   ! Nitrogênio da serapilheira
      real(r_4),dimension(pl) :: pl_phosphorus = 0.0 ! FORFI do serapilheira
      real(r_4),dimension(ps) :: ps_nitrogen = 0.0   ! & so forth
      real(r_4),dimension(ps) :: ps_phosphorus = 0.0

      real(r_4) :: leaf_n  ! Mass of Nutrients in MO comming from Vegetation  g(Nutrient m⁻²)
      real(r_4) :: froot_n
      real(r_4) :: wood_n
      real(r_4) :: leaf_p
      real(r_4) :: froot_p
      real(r_4) :: wood_p

      real(r_4) :: water_modifier ! Multiplicator for water influence on C decay

      real(r_4) :: frac1,frac2

      real(r_4),dimension(pl+ps) :: het_resp, cdec
      real(r_4),dimension(pl+ps) :: aux_ratio_n, aux_ratio_p
      real(r_4),dimension(pl+ps) :: nutri_min_n, nutri_min_p


      ! Need to define partition of C, N & P among Soil Pools (A function based on Lignin content c(w)ould be great)
      frac1 = 0.8
      frac2 = 1.0 - frac1

      snr = 0.0

      ! Turnover Rates  == residence_time⁻¹ (years⁻¹)
      ! Change it in future (Parametrize from literarure)
      tr_c(1) = 10           ! litter I   (1)
      tr_c(2) = 30.0        ! litter II  (2)
      tr_c(3) = 300.0       ! soil   I   (3)
      tr_c(4) = 700.0      ! soil   II  (4)

      ! find nutrient mass/area) : litter fluxes[ML⁻²] * litter nutrient ratios
      ! (lnr) [MM⁻¹]
      leaf_n = leaf_l * lnr(1) ! g(nutrient) m-2
      froot_n = root_l * lnr(2)
      wood_n = cwd * lnr(3)
      leaf_p = leaf_l * lnr(4)
      froot_p = root_l * lnr(5)
      wood_p = cwd * lnr(6)

      ! FIRST OF ALL calculate dacay from pools
      ! CARBON DECAY
      water_modifier = water_effect(water_sat)
      do index = 1,4
         if(index .lt. 3) then
            ! FOR THE 2 LITTER POOLS
            cdec(index) = scarbon_decayment(q10,tsoil,cl(index),tr_c(index)) * water_modifier
         else
            ! FOR THE 3 CARBON POOLS
            cdec(index) = scarbon_decayment(q10,tsoil,cs(index-2),tr_c(index)) * water_modifier
         endif
      enddo

      ! partitioning material coming from vegetation
      ! filling C (litter I-II; soil I) pools with incoming material
      !LITTER I
      aux1 = (frac1 * leaf_l) + (frac1 * root_l)                             ! INcoming Carbon
      aux2 = cdec(1) * clit_atm                                              ! Carbon lost to ATM
      aux3 = cdec(1) - aux2                                                  ! Carbon going to cl_out(2)
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
            aux_ratio_n(index) = nmass_org(index)/cl_out(index) ! g(N)g(C)-1
            aux_ratio_p(index) = pmass_org(index)/cl_out(index) ! g(P)g(C)-1
         else
            aux_ratio_n(index) = nmass_org(index)/cs_out(index-2) ! g(N)g(C)-1
            aux_ratio_p(index) = pmass_org(index)/cs_out(index-2) ! g(P)g(C)-1
         end if
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


      ! USE NUTRIENT RATIOS AND HET_RESP TO CALCULATE MINERALIZED NUTRIENTS
      do index=1,4
         nutri_min_n(index) = het_resp(index)*aux_ratio_n(index)
         nutri_min_p(index) = het_resp(index)*aux_ratio_p(index)
      enddo

      ! UPDATE N and P in SOIL POOLS
      do index=1,4
         nmass_org(index) = nmass_org(index) - nutri_min_n(index)
         pmass_org(index) = pmass_org(index) - nutri_min_p(index)
      enddo

      hr = sum(het_resp)


      ! FEATURES TO BE IMPLEMENTED

      ! P weathering +
      ! P biochemical mineralization +
      ! P release "de-sorption" +
      ! P immobilization
      ! P leaching -
      ! P sorption I & II -
      ! P occlusion -

      ! N deposition +
      ! N fixation +
      ! N immobilization -
      ! N leaching -
      ! N degassing-volatilization -

!       mineral_n_glob = mineral_n_glob + (sum(nutri_min_n,&
!            &mask=.not. isnan(nutri_min_n)) * 1e-3) ! - kg m⁻²
!       labile_p_glob = labile_p_glob + (sum(nutri_min_p,&
!            &mask=.not. isnan(nutri_min_p)) * 1e-3) ! - kg m⁻²

      ! print *, 'c AUX', carbon_aux
      ! print *, 'n mineral', nutri_min_n
      ! print *, 'p labila', nutri_min_p

      ! print *, 'carbon3 _pools'
      ! print *, 'labile_p_glob ---.>',labile_p_glob
      ! print *, 'mineral_n_glob ---.>',mineral_n_glob
      ! print *, 'hr ---.>',hr

   end subroutine carbon3

end module soil_dec
