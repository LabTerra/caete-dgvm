program test_carbon3
   use types
   use utils
   use global_par
   use soil_dec

   implicit none

   print *, "Testing carbon decayment funtion"
   call test_scarbon_t()
   
   print *,
   print *,
   print *, "Testing water influence funtion"
   call test_water_function()
   
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
         result(j) = scarbon_decayment(q101, temp_range(j), c, res_time)
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

end program test_carbon3 