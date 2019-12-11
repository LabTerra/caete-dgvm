program test_carbon3
   use types
   use utils
   use global_par
   use photo
   use soil_dec

   implicit none

   ! print *, "Testing carbon decayment funtion"
   ! call test_scarbon_t()

   ! print *,
   ! print *,
   ! print *, "Testing water influence funtion"
   ! call test_water_function()


   ! print *,
   ! print *,
   ! print *, "Testing/debugging CARBON3"

   ! call test_c3()

   print *,
   print *,
   print *, "Testing/debugging Allocation"

   call test_alloc()


   contains

   ! subroutine template(arg1,  arg2)
   !    implicit none
   !    real(r_4),intent(in) :: arg1
   !    real(r_4),intent(out) ::  arg2
   !    arg2 = arg1 - 1.0
   ! end subroutine template

   ! function template1(arg) result(retval)
   !    implicit none
   !    real(r_4) :: arg
   !    real(r_4) :: retval

   ! end function template1

   subroutine test_scarbon_t()
      ! TEST FOR carbon_dacayment function
      ! DEfinition: arguments : (q10,tsoil,c,residence_time) result(decay)
      ! TYPES
      ! real(r_4),intent(in) :: q10              ! constant ~1.4
      ! real(r_4),intent(in) :: tsoil            ! Soil temperature °C
      ! real(r_4),intent(in) :: c                ! Carbon content per area g(C)m-2
      ! real(r_4),intent(in) :: residence_time   ! Pool turnover rate
      ! real(r_4) :: decay ! ML⁻²

      integer(i_4) :: i, j

      real(kind=r_4) :: q101 = q10,   &
                     & c = 20.0,    &
                     & res_time = 20

      real(kind=r_4), allocatable, dimension(:) :: temp_range, result

      allocate(temp_range(50))
      allocate(result(50))

      call linspace(from=-5.0, to=45.0, array=temp_range)

      do j = 1,50
         result(j) = carbon_decay(q101, temp_range(j), c, res_time)
      end do

      do i = 1, 50
         print *, result(i)
      end do

      deallocate(temp_range)
      deallocate(result)

   end subroutine test_scarbon_t

   !---------------------------------------------------------------------
   ! TEST FOR water influence on Carbon decaiment
   subroutine test_water_function()
      implicit none
      integer(i_4) :: i
      real(r_4),dimension(50) :: water_range

      call linspace(from=0.0, to=1.0, array=water_range)

      do i = 1,50
         print*, water_effect(water_range(i))
      enddo

   end subroutine test_water_function


   !---------------------------------------------------------------------
   ! TEST CARBON3

   subroutine test_c3()

      integer(i_4) :: index, j
      real(r_4) :: soilt=23.0, water_s=0.8, ll=0.0000001, lf=0.0000001, lw=0.0000001
      real(r_4), dimension(6) :: lnr = (/0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.00001/)
      real(r_4), dimension(2) :: cl = 0.0, cs = 0.0, cl_out = 0.0, cs_out = 0.0
      real(r_4), dimension(8) :: snr = 0.0
      real(r_4) :: hr, nupt, pupt

      pupt = 0.05
      nupt = 0.02

      do index = 1,1000000
         call carbon3(soilt,water_s, ll, lw, lf, lnr, cl, cs, nupt, pupt, cl_out, cs_out, snr, hr)
         do j = 1,2
            cs(j) = cs_out(j)
            cl(j) = cl_out(j)
         end do
      end do

      print *, snr,"<- snr"
      print *, hr,"<- hr"
      print *, cl,"<- cl"
      print *, cs,"<- cs"

   end subroutine test_c3

   subroutine test_alloc()

      ! INPUTS
      real(r_4),dimension(ntraits) :: dt = (/12.799142793925565, 0.0,4.532783508300781,44.61568069458008,&
                                           &  6.493612289428711, 0.5559999999999985,0.21599999999999978,&
                                           &  0.22799999999999973, 0.0, 0.05658099806223465,0.0749577527005275,&
                                           &  0.02557946262317787, 0.007156888889925779,0.007339313026394539,&
                                           &  0.00636147002287577/)  ! PLS attributes
      real(r_4) :: npp = 0.1  ! npp (KgC/m2/yr) from assimilation process
      real(r_8) :: scl1 = 0.5d0 ! previous day carbon content on leaf compartment (KgC/m2)
      real(r_8) :: sca1 = 7.0d0 ! previous day carbon content on aboveground woody biomass compartment(KgC/m2)
      real(r_8) :: scf1 = 0.5d0! previous day carbon content on fine roots compartment (KgC/m2)
      real(r_8) :: nmin = 0.00002d0 ! N in mineral N pool(kg m-2)
      real(r_8) :: plab  = 0.00001d0 ! P in labile pool (kg m-2)
      real(r_8),dimension(3) :: storage = (/0.0d0, 0.0d0, 0.0d0/)! Three element array- storage pool([C,N,P]) g m-2

      ! OUTPUTS
      real(r_8),dimension(3) :: storage_out
      real(r_8) :: scl2 ! final carbon content on leaf compartment (KgC/m2)
      real(r_8) :: sca2 ! final carbon content on aboveground woody biomass compartment (KgC/m2)
      real(r_8) :: scf2 ! final carbon content on fine roots compartment (KgC/m2)
      real(r_8) :: cwd  ! coarse wood debris (to litter)(C) g m-2
      real(r_8) :: root_litter ! to litter g(C) m-2
      real(r_8) :: leaf_litter ! to litter g(C) m-2
      real(r_8) :: nuptk ! N plant uptake g(N) m-2
      real(r_8) :: puptk ! P plant uptake g(P) m-2
      real(r_8),dimension(6) :: litter_nutrient_ratios ! [(lln2c),(rln2c),(cwdn2c),(llp2c),(rlp2c),(cwdp2c)]
      logical(l_1) :: end_pls_day ! ABORT MISSION SIGN

      integer(l_1) :: index


      do index = 1,100

         call allocation(dt, npp, nmin,plab,scl1,sca1,scf1,storage,&
         &storage_out,scl2,sca2,scf2,leaf_litter,cwd,root_litter,&
         &nuptk,puptk,litter_nutrient_ratios,end_pls_day)

         scl1 = scl2
         sca1 = sca2
         scf1 = scf2
         storage = storage_out

         print *, scl1, sca1, scf1, "-> C pools"
         print *, storage, "      -> Storage"
         print *, puptk, "        -> PUPTK"
         print *, nuptk, "        -> NUPTK"
      end do

   end subroutine test_alloc

end program test_carbon3
