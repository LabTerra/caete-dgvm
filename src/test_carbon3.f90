program test_carbon3
   use types
   use global_par
   use soil_dec

   implicit none

   
   ! TEST FOR carbon_dacayment function
   ! DEfinition: arguments : (q10,tsoil,c,residence_time) result(decay)
   ! TYPES
   ! real(r_4),intent(in) :: q10              ! constant ~1.4
   ! real(r_4),intent(in) :: tsoil            ! Soil temperature °C
   ! real(r_4),intent(in) :: c                ! Carbon content per area g(C)m-2
   ! real(r_4),intent(in) :: residence_time   ! Pool turnover rate
   ! real(r_4) :: decay ! ML⁻²
   
   integer(i_4) :: i, j
   
   real(kind=r_4) :: q10 = 2.1,   &
                   & tsoil = 23.0,&
                   & c = 20.0,    &
                   & res_time = 20

   real(kind=r_4) :: decayment
   !real(kind= r_4), external :: scarbon_decayment 
   real(kind=r_4), allocatable, dimension(:) :: temp_range, result
   
   allocate(temp_range(50))
   allocate(result(50))
   
   do i = -5,45
      temp_range(i + 5) = i 
   end do
  
   do j = 1,50
      result(j) = scarbon_decayment(q10, temp_range(j), c, res_time)
   end do

   do i = 1, 50
      print *, result(i)
   end do

end program test_carbon3
