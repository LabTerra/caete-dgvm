! Author: João Paulo Darela FIlho- joao.darela@unesp.br


module soil_dec
   ! This module represents a global struct (This will only work if isolated in a sinlge process)
   ! To be used only in caete_dyn  and daily_budget subroutines
   use types
   use global_par
   
   implicit none
   private
   
   ! ## STATE VARIABLES
   real(r_4),public,dimension(4) :: carbon_aux = 0.0    ! Variável Auxiliar
   real(r_4),public,dimension(2) :: litter_carbon = 0.0 ! Litter carbon pools
   real(r_4),public,dimension(3) :: soil_carbon = 0.0   ! soil carbon pools

   real(r_4),public,dimension(4) :: soil_nr    ! Soil pools Nutrient Ratios (2N + 2P)
   real(r_4),public,dimension(4) :: litt_nr    ! litter pools Nutrient Ratios (2N + 2P)
!
!    real(r_4),public,dimension(2) :: litter_p
!    real(r_4),public,dimension(2) :: litter_n
!    real(r_4),public,dimension(3) :: soil_p
!    real(r_4),public,dimension(3) :: soil_n

   ! These are global variables that are initialized in caete_init.
   ! New inputs -- Estimated from literature for Manaus Region
   real(r_8),public :: mineral_n_glob = 3.775999d-4 ! kg m-2 Xu et al. 2013 ?
   real(r_8),public :: labile_p_glob = 2.4299955d-4  ! kg m-2 Yang et al., 2013

   !  real(r_4),public :: occluded_p
   !  real(r_4),public :: sorbed_p
   !  real(r_4),public :: ssorbed_p
   logical(l_1),public :: last_run = .false.

   ! These are global variables that are initialized in caete_init.
   public :: process_id, carbon3
   public :: scarbon_decayment, water_effect
   public :: set_var


contains


   subroutine set_var(arg1, arg2)

      real(kind=r_8),intent(   in) :: arg1
      real(kind=r_8),intent(inout) :: arg2

      arg2 = arg1

   end subroutine set_var


   ! Identify process number
   function process_id() result(ipid)
      integer(r_4) :: ipid

      ipid = getpid()

   end function process_id


   !Based on carbon decay implemented in JeDi and JSBACH - Pavlick et al. 2012
   function scarbon_decayment(q10_in,tsoil,c,residence_time) result(decay)

      real(r_4),intent(in) :: q10_in              ! constant ~1.4
      real(r_4),intent(in) :: tsoil            ! Soil temperature °C
      real(r_4),intent(in) :: c                ! Carbon content per area g(C)m-2
      real(r_4),intent(in) :: residence_time   ! Pool turnover rate
      real(r_4) :: decay ! ML⁻²

      if(c .le. 0.0) then
         decay = 0.0
         return
      endif

      decay = (q10_in**((tsoil - 20.0) / 10.0)) * (c/residence_time)

   end function scarbon_decayment
   


   function water_effect(theta) result(retval)
      ! Implement the Moyano function based on soil water content 
      real(r_4),intent(in) :: theta  ! Volumetric soil water content (cm³ cm⁻³)
      real(r_4),parameter :: k_a = 3.11, k_b = 2.42 
      real(r_4) :: retval

      retval= (k_a * theta) - (k_b * theta**2)
   
   end function water_effect


   subroutine carbon3(tsoil,water_sat,leaf_l, cwd, root_l, lnr, cl, cs, cl_out, cs_out, soil_nr, hr)
      ! this one wastes 132 bits of primary memory per process

      integer(i_4),parameter :: pl=2,ps=3
      integer(i_4) :: index

      !     Variables
      !     =========
      !     Inputs
      !     ------
      real(r_4),intent(in) :: tsoil, water_sat       ! soil temperature (oC); soil water relative content (dimensionless)
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
      real(r_4),dimension(4), intent(out) :: soil_nr ! Soil pools Nutrient to C ratios 
      real(r_4),intent(out) :: hr                    !Heterotrophic (microbial) respiration (gC/m2/day)

      !TODO ! Insert output: Total mineralized N and P
      !     Internal
      real(r_4),dimension(pl+ps) :: tr_c
      real(r_4),dimension(pl) :: pl_nitrogen = 0.0   ! Nitrogênio da serapilheira
      real(r_4),dimension(pl) :: pl_phosphorus = 0.0 ! FORFI do serapilheira
      real(r_4),dimension(ps) :: ps_nitrogen = 0.0   ! & so forth
      real(r_4),dimension(ps) :: ps_phosphorus = 0.0

      real(r_4) :: leaf_n2c  ! Mass of Nutrients in MO comming from Vegetation  g(Nutrient m⁻²)
      real(r_4) :: froot_n2c
      real(r_4) :: wood_n2c
      real(r_4) :: leaf_p2c
      real(r_4) :: froot_p2c
      real(r_4) :: wood_p2c
      real(r_4) :: water_modifier
      real(r_4) :: frac1,frac2
      real(r_4),dimension(pl+ps) :: het_resp, cdec
      real(r_4),dimension(pl+ps) :: aux_ratio_n, aux_ratio_p
      real(r_4),dimension(pl+ps) :: nutri_min_n, nutri_min_p

      ! PARAMETERS
      ! Fraction of C that is broken by soil microbial activity and goes to Atmosphere (A.K.A. Heterothrophic respiration )
      ! ! From Pavlick et al. 2012 - JeDi/ Raddatz et al. 2007 - JSBACH
      real(r_4),parameter :: clit_atm = 0.77
      real(r_4),parameter :: cwd_atm = 0.2

      ! Need to define partition of C, N & P among Soil Pools (A function based on Lignin content c(w)ould be great)
      frac1 = 0.7
      frac2 = 1.0 - frac1


      soil_nr = 0.0
      ! Turnover Rates  == residence_time⁻¹ (years⁻¹)
      ! Change it in future (Parametrize from literarure)
      tr_c(1) = 5    ! litter I   (1)
      tr_c(2) = 30.0   ! litter II  (2)
      tr_c(3) = 800.0   ! soil   I   (3)
      tr_c(4) = 7000.0  ! soil   II  (4)
      tr_c(5) = 0.0 ! soil   III (5) 4B

      ! find nutrient mass/area) : litter fluxes[ML⁻²] * litter nutrient ratios
      ! (lnr) [MM⁻¹]
      leaf_n2c = leaf_l * lnr(1) ! g(nutrient) m-2
      froot_n2c = root_l * lnr(2)
      wood_n2c = cwd * lnr(3)
      leaf_p2c = leaf_l * lnr(4)
      froot_p2c = root_l * lnr(5)
      wood_p2c = cwd * lnr(6)

      ! FIRST OF ALL calculate dacay from pools
      ! CARBON DECAY
      water_modifier = water_effect(water_sat)
      do index = 1,4
         if(index .lt. 3) then
            ! FOR THE 2 LITTER POOLS
            cdec(index) = scarbon_decayment(q10,tsoil,cl(index),tr_c(index)) * water_modifier
         else
            ! FOR THE 3 CARBON POOLS
            cdec(index) = scarbon_decayment(q10,tsoil,cs(index-2),tr_c(index))  * water_modifier
         endif
      enddo
      cdec(5) = 0.0 ! + 4B

      ! partitioning material coming from vegetation
      ! filling C (litter I-II; soil I) pools with incoming material
      !LITTER I & II

      ! CONTROL HERE THE CARBON FLOW COMING FROM BAD PLSs (DURING SPINUP)
      cl_out(1) = (cl(1) - cdec(1)) + (frac1 * leaf_l) + (frac1 * root_l)
      cl_out(2) = (cl(2) - cdec(2)) + (frac2 * leaf_l) + (frac2 * root_l)
      cl_out(2) = cl_out(2) + (cwd * frac2)
      !SOIL Ie
      cs_out(1) = (cs(1) - cdec(3)) + (frac1 * cwd)
      cs_out(2) = ((1 - clit_atm) * (cdec(1) + cdec(2))) + &
           & ((1 - cwd_atm) * cdec(3)) - cdec(4)
      cs_out(3) = 0.0 ! cs_out(3) is empty ! + 4B
      ! print *, 'cl after partitioning->', cl
      ! print *, 'cs after partitioning->', cs

      ! Calculate the nutrint contents in each litterI-II/soilI pool
      pl_nitrogen(1)   = (leaf_n2c * frac1)+(froot_n2c * frac1)  ! g(N)m-2
      pl_nitrogen(2)   = (leaf_n2c * frac2)+(froot_n2c * frac2)+(wood_n2c *frac2)
      ps_nitrogen(1)   = wood_n2c * frac1
      pl_phosphorus(1) = (leaf_p2c * frac1)+(froot_p2c * frac1)
      pl_phosphorus(2) = (leaf_p2c * frac2)+(froot_p2c * frac2)+(wood_p2c*frac2)
      ps_phosphorus(1) = wood_p2c * frac1

      ! System Memory !WARNING! carbon aux is a global variable
      ! Storage N and P that are bonded with C in the SLOW (do not receives MO from vegetation)
      ! SLOW POOL only receives MO(recalcitrant) from L1, L2 and S1 pools
      ps_nitrogen(2) = carbon_aux(1)
      ps_nitrogen(3) = 0.0 ! + 4B ! O QUINTO ELEMENTO - UNUNSED VARIABLE
      ps_phosphorus(2) = carbon_aux(3)
      ps_phosphorus(3) = 0.0 ! + 4B ! O QUINTO ELEMENTO - UNUNSED VARIABLE

      !HRESP
      het_resp(1) = clit_atm * cdec(1)
      het_resp(2) = clit_atm * cdec(2)
      het_resp(3) = cwd_atm * cdec(3)
      het_resp(4) = cdec(4)
      het_resp(5) = 0.0 ! + 4B

      ! THIS SECTION - Mineralized nutrients =====================================================
      ! The amounts of Minerilized nutrients are dependent on N:C and P:C mass ratios of soil pools
      ! NUTRIENT RATIOS in SOIL
      do index=1,2
         ! FIRST TO LITTER POOLS (2x)
         if(cl_out(index) .gt. 0.0) then
            ! Nutrient:Carbon Mass Ratio of litter pools
            aux_ratio_n(index) = pl_nitrogen(index)/cl_out(index)   ! g(N)g(C)-1
            aux_ratio_p(index) = pl_phosphorus(index)/cl_out(index) ! g(P)g(C)-1
         else
            aux_ratio_n(index) = 0.0
            aux_ratio_p(index) = 0.0
         endif
      enddo
      ! **** To GLOBAL VARS
      litt_nr(1) = aux_ratio_n(1)
      litt_nr(2) = aux_ratio_n(2)
      litt_nr(3) = aux_ratio_p(1)
      litt_nr(4) = aux_ratio_p(2)

      do index=3,4
         ! THEN TO SOIL POOLS (2X)
         if(cs_out(index-2) .gt. 0.0) then
            ! Nutrient:Carbon Mass Ratio of soil pools
            aux_ratio_n(index) = ps_nitrogen(index-2)/cs_out(index-2)   ! g(N)g(C)-1
            aux_ratio_p(index) = ps_phosphorus(index-2)/cs_out(index-2) ! g(P)g(C)-1
         else
            aux_ratio_n(index) = 0.0
            aux_ratio_p(index) = 0.0
         endif
      enddo
      ! ==== To GLOBAL VARS
      soil_nr(1) = aux_ratio_n(3)
      soil_nr(2) = aux_ratio_n(4)
      soil_nr(3) = aux_ratio_p(3)
      soil_nr(4) = aux_ratio_p(4)

      aux_ratio_p(5) = 0.0 ! + 4B

      ! USE NUTRIENT RATIOS AND HET_RESP TO CALCULATE MINERALIZED NUTRIENTS
      do index=1,4
         nutri_min_n(index) = het_resp(index)*aux_ratio_n(index)
         nutri_min_p(index) = het_resp(index)*aux_ratio_p(index)
      enddo
      nutri_min_n(5) = 0.0 ! + 4B
      nutri_min_p(5) = 0.0 ! + 4B

      ! UPDATE N and P in SOIL POOLS
      do index=1,2
         ! LITTER
         pl_nitrogen(index) = pl_nitrogen(index) - nutri_min_n(index)
         pl_phosphorus(index) = pl_phosphorus(index) - nutri_min_p(index)
      enddo

      do index=1,2
         ! SOIL
         ps_nitrogen(index) = ps_nitrogen(index) - nutri_min_n(index + 2)
         ps_phosphorus(index) = ps_phosphorus(index) - nutri_min_p(index + 2)
      enddo

      ! DO SOMETHING WITH THESE UNUNSED VARIABLES SO THE COMPILER DONT WARNS YOU
      ps_nitrogen(3) = 0.0
      ps_phosphorus(3) = 0.0

      ! calculate final state
      ! SOIL II and III nutrient pools
      carbon_aux(1) = ps_nitrogen(2)
      carbon_aux(2) = 0.0
      carbon_aux(3) = ps_phosphorus(2)
      carbon_aux(4) = 0.0

      if(last_run) then
         continue
         ! fill soil_dec_variables
      endif

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

      hr = sum(het_resp)

      ! print *, 'c AUX', carbon_aux
      ! print *, 'n mineral', nutri_min_n
      ! print *, 'p labila', nutri_min_p

      ! print *, 'carbon3 _pools'
      ! print *, 'labile_p_glob ---.>',labile_p_glob
      ! print *, 'mineral_n_glob ---.>',mineral_n_glob
      ! print *, 'hr ---.>',hr

   end subroutine carbon3

end module soil_dec