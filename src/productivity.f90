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

! contacts :: David Montenegro Lapola <lapoladm ( at ) gmail.com> 


module productivity
  implicit none
  private

  public :: prod


contains

  subroutine prod(dt,light_limit,temp,p0,w,ipar,rh,emax,cl1_prod,&
       & ca1_prod,cf1_prod,beta_leaf,beta_awood,beta_froot,sto1,ph,ar,&
       & nppa,laia,f5,vpd,rm,rg,rc,wue,c_defcit,vm_out,sla,sto2)

    use types
    use global_par
    use photo
    use water

!Input
!-----     
    real(r_4),dimension(ntraits),intent(in) :: dt ! PLS data 
    real(r_4), intent(in) :: temp                 !Mean monthly temperature (oC)
    real(r_4), intent(in) :: p0                   !Mean surface pressure (hPa)
    real(r_4), intent(in) :: w                    !Soil moisture kg m-2
    real(r_4), intent(in) :: ipar                 !Incident photosynthetic active radiation (w/m2)
    real(r_4), intent(in) :: rh,emax !Relative humidity/MAXIMUM EVAPOTRANSPIRATION
    real(r_8), intent(in) :: cl1_prod, cf1_prod, ca1_prod        !Carbon in plant tissues (kg/m2)
    real(r_8), intent(in) :: beta_leaf            !npp allocation to carbon pools (kg/m2/day)
    real(r_8), intent(in) :: beta_awood
    real(r_8), intent(in) :: beta_froot
    real(r_8), dimension(3), intent(in) :: sto1
    logical(l_1), intent(in) :: light_limit                !True for no ligth limitation

!     Output
!     ------
    real(r_4), intent(out) :: ph                   !Canopy gross photosynthesis (kgC/m2/yr)
    real(r_4), intent(out) :: rc                   !Stomatal resistence (not scaled to canopy!) (s/m)
    real(r_4), intent(out) :: laia                 !Autotrophic respiration (kgC/m2/yr)
    real(r_4), intent(out) :: ar                   !Leaf area index (m2 leaf/m2 area)
    real(r_4), intent(out) :: nppa                 !Net primary productivity (kgC/m2/yr) 
    real(r_4), intent(out) :: vpd            
    real(r_4), intent(out) :: f5                   !Water stress response modifier (unitless) 
    real(r_4), intent(out) :: rm                   !autothrophic respiration (kgC/m2/day)
    real(r_4), intent(out) :: rg 
    real(r_4), intent(out) :: wue
    real(r_4), intent(out) :: c_defcit     ! Carbon deficit gm-2 if it is positive, aresp was greater than npp + sto2(1) 
    real(r_8), intent(out) :: sla          !specific leaf area (m2/kg)
    real(r_8), intent(out) :: vm_out
    real(r_8), dimension(3), intent(out) :: sto2
!     Internal
!     --------

    real(r_4) :: tleaf,awood            !leaf/wood turnover time (yr)
    real(r_4) :: g1
    real(r_4) :: c4

    real(r_4) :: n2cl
    real(r_4) :: n2cl_resp
    real(r_4) :: n2cw_resp
    real(r_4) :: n2cf_resp
    real(r_4) :: p2cl
    integer(i_4) :: c4_int

    real(r_4) :: f1       !Leaf level gross photosynthesis (molCO2/m2/s)
    real(r_4) :: f1a      !auxiliar_f1
    logical :: no_cell = .false.

!getting pls parameters

    if(((cl1_prod .lt. cmin) .and. (cf1_prod .lt. cmin))) then !  
! Then PLS 'Die'
       no_cell = .true.
       goto 999
    endif

    g1  = dt(1)
    tleaf = dt(3)
    awood = dt(7)
    c4  = dt(9)
    n2cl = dt(10)
    n2cl_resp = n2cl    
    n2cw_resp = dt(11)
    n2cf_resp = dt(12)
    p2cl = dt(13)

! open (unit=12324,file='dt_dentro.txt',status='unknown',&
!         &form='formatted',position='append', action='write')

! write (12324,*) dt

! close(12324)

    n2cl = real(n2cl * (cl1_prod * 1e3), r_4) ! N in leaf g m-2
    p2cl = real(p2cl * (cl1_prod * 1e3), r_4) ! P in leaf g m-2
    c4_int = nint(c4) 

    if(debug) then    
       write(1234,*) '-----Message from productivity-----------------'
       write(1234,*) '-----------------------------------------------'
    endif

!     ==============
!     Photosynthesis 
!     ==============
! rate (molCO2/m2/s)
! subroutine photosynthesis_rate(temp,p0,ipar,ll,c3,nbio,pbio,&
!                           & leaf_turnover,storage_1,f1ab,vm,storage_2)
! print *, n2cl, "N gm-2"
! print *, p2cl, "P gm-2"

    call photosynthesis_rate(temp,p0,ipar,light_limit,c4_int,n2cl,&
         & p2cl,tleaf,sto1,f1a,vm_out,sto2)


    if(debug) then 
       write(1234,*) 'f1a -->',f1a
       if(f1a < 0.0) then
          print *, 'f1a less than 0 aborting'
          call abort
       endif
    endif

    if(debug) then
       write(1234,*)
    endif

! VPD
!========
    vpd = vapor_p_defcit(temp,rh)

    if(debug) then
       write(1234,*) 'VPD -->', vpd
    endif

!Stomatal resistence
!===================
    rc = canopy_resistence(vpd, f1a, g1)

! Novo calculo da WUE

    wue = water_ue(f1a, rc, p0, vpd)

!     Water stress response modifier (dimensionless)
!     ----------------------------------------------
    f5 =  water_stress_modifier(w, cf1_prod, rc, emax)

    if(debug) then
       write(1234,*) 'f5c ->', f5
    endif

!     Photosysthesis minimum and maximum temperature
!     ----------------------------------------------

    if ((temp.ge.-10.0).and.(temp.le.50.0)) then
       f1 = f1a * f5 ! :water stress factor ! Ancient floating-point underflow spring (from CPTEC-PVM2)
    else
       f1 = 0.0      !Temperature above/below photosynthesis windown
    endif

    if(debug) then
       write(1234,*) 'f1 (after f5)->', f1
    endif


!     Leaf area index (m2/m2)
    sla = spec_leaf_area(tleaf)
    laia = real(leaf_area_index(cl1_prod, sla), r_4)
! laia = real(f_four(90,cl1_prod,sla), r_4)          ! sunlai
! laia = laia + (real(f_four(20,cl1_prod,sla), r_4)) ! shadelai


    if(debug) then
       write(1234,*) 'tleaf ->', tleaf
       write(1234,*) 'SLA ->', sla
       write(1234,*) 'LAI ->', laia
       write(1234,*) 'betaleaf ->', beta_leaf
       write(1234,*) 'betawood ->', beta_awood
       write(1234,*) 'betaroot ->', beta_froot
    endif

!     Canopy gross photosynthesis (kgC/m2/yr)
!     =======================================x
    ph =  gross_ph(f1,cl1_prod,sla)       ! kg m-2 year-1


    if(debug) then
       write(1234,*) 'ph ->', ph
    endif
!     Autothrophic respiration
!     ========================
!     Maintenance respiration (kgC/m2/yr) (based in Ryan 1991)
    rm = m_resp(temp,cl1_prod,cf1_prod,ca1_prod,sto2&
         &,n2cl_resp,n2cw_resp,n2cf_resp,awood)

! c     Growth respiration (KgC/m2/yr)(based in Ryan 1991; Sitch et al.
! c     2003; Levis et al. 2004)         
    rg = g_resp(beta_leaf,beta_awood, beta_froot,awood) 

    if (rg.lt.0) then
       rg = 0.0
    endif

!     c Autotrophic (plant) respiration -ar- (kgC/m2/yr)
!     Respiration minimum and maximum temperature
!     -------------------------------------------     
    if ((temp.ge.-10.0).and.(temp.le.50.0)) then
       ar = rm + rg
    else
       ar = 0.0               !Temperature above/below respiration windown
    endif

!     Net primary productivity(kgC/m2/yr)
!     ====================================
    nppa = ph - ar
! this operation affects the model mass balance
! If ar is bigger than ph, what is the source or respired C?

    if(ar .gt. ph) then
       c_defcit = ((ar - ph) * 2.73791) ! tranform kg m-2 year-1 in  g m-2
       nppa = 0.0
       if(c_defcit .gt. real(sto2(1),r_4))then
          c_defcit =  c_defcit - real(sto2(1),r_4)
          sto2(1) = 0.0
       else ! c_defcit  is .le. sto2(1)
          c_defcit = 0.0
          sto2(1) = sto2(1) - c_defcit
       endif
    else
       c_defcit = 0.0 
    endif

    no_cell = .false.


    if(debug) then
       write(1234,*) 'sto2 ->', sto2
       write(1234,*) 'nppa ->', nppa
       write(1234,*) 'ar ->', ar
       write(1234,*) 'cdef ->', c_defcit
    endif


    if(debug) then
       write(1234,*) '-----END Message from productivity--- ---------'
       write(1234,*) '-----------------------------------------------'
    endif

999 continue
    if(no_cell) then
       ph = 0.0                   
       rc = 0.0                   
       laia = 0.0                 
       ar = 0.0                   
       nppa = 0.0                  
       vpd = 0.0            
       f5 = 0.0                   
       rm = 0.0             
       rg = 0.0
       wue = 0.0
       c_defcit = 0.0
       sla = 0.0          
       vm_out = 0.0
       sto2 = 0.0
    endif

  end subroutine prod

end module productivity

